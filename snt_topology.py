#!/usr/bin/env python3
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import pathlib, glob
from matplotlib import cm
from matplotlib.colors import ListedColormap
from sklearn.linear_model import LinearRegression
from scipy.stats import kstest
from scipy.stats import pearsonr
from scipy.stats import chisquare
from scipy.interpolate import splprep, splev
from scipy.integrate import simpson
from typing import TypeVar
from typing import Tuple
from typing import Union
Self = TypeVar('Self', bound = 'SNT_graph')

class SNT_graph:
    def __init__(self: Self,
                 df: pd.DataFrame,
                 xy_res: float = .065,
                 z_res: float = .2,
                 min_len: float = 1.5,
                 min_cycle: float = 5.,
                 ) -> Self:
        self.xy_res = xy_res
        self.z_res = z_res
        self.min_len = min_len
        self.min_cycle = min_cycle
        self.G = nx.Graph()
        self.add_edges_from_df(df)
        self.merge_junctions()
        self.graph_statistics()

    def add_edges_from_df(self: Self,
                          df: pd.DataFrame,
                          init_path: int = 2,
                          max_n_paths: int = 100
                          ) -> None:
        ix = np.where(df['Name'] == 'Apex')[0]
        self.apex_zyx = df[['Z','Y','X']].iloc[ix].values[0]
        scaling = np.array([self.z_res, self.xy_res, self.xy_res])
        for i in range(init_path, max_n_paths):
            ix = np.where(df['Name'] == f'Path ({i})-XY-0001')[0]
            if ix.size > 0:
                for ii in range(ix.size - 1):
                    ix0 = df.loc[ix[ii], ['Z','Y','X']].values
                    ix1 = df.loc[ix[ii + 1], ['Z','Y','X']].values
                    w01 = np.linalg.norm((ix1 - ix0) * scaling)
                    self.G.add_edge(str(ix0), str(ix1), weight = w01, pll_branch = (i == 2))
        self.H = self.G.copy()
        self.simplify_and_prune()

    def merge_junctions(self: Self) -> None:
        self.remove_single_gaps()
        self.remove_double_gaps()
        self.remove_short_cycles()
        self.remove_triple_gaps()
        self.simplify_and_prune()

    def remove_single_gaps(self: Self) -> None:
        # Merge junctions across zero-node gaps
        for node in list(self.G.nodes()):
            if self.G.degree(node) > 2:
                nb = list(self.G.neighbors(node))
                for nb1 in nb:
                    if self.G.has_node(node):
                        if self.G.degree(nb1) > 2:
                            for nb2 in nb:
                                w0 = self.G.edges[node, nb1]['weight']
                                if nb2 is not nb1:
                                    w0 += self.G.edges[node, nb2]['weight']
                                    b = self.G.edges[node, nb2]['pll_branch']
                                    self.G.add_edge(nb1, nb2, weight = w0, pll_branch = b)
                            self.G.remove_node(node)

    def remove_double_gaps(self: Self) -> None:
        # Merge junctions across one-node gaps
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 2:
                nb = list(self.G.neighbors(node))
                nb_deg = list(self.G.degree(nb))
                if nb_deg[0][1] > 2 and nb_deg[1][1] > 2:
                    for nb1 in list(self.G.neighbors(nb[0])):
                        w0 = self.G.edges[node, nb[0]]['weight']
                        w0 += self.G.edges[node, nb[1]]['weight']
                        if nb1 is not node:
                            w0 += self.G.edges[nb[0], nb1]['weight']
                            b = self.G.edges[nb[0], nb1]['pll_branch']
                            self.G.add_edge(nb[1], nb1, weight = w0, pll_branch = b)
                    self.G.remove_nodes_from([node, nb[0]])

    def remove_triple_gaps(self: Self) -> None:
        # Merge junctions across two-node gaps
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 2:
                nb = list(self.G.neighbors(node))
                nb1, nb0 = None, None
                if self.G.degree(nb[0]) == 2:
                    nb1, nb0 = nb[0], nb[1]
                elif self.G.degree(nb[1]) == 2:
                    nb1, nb0 = nb[1], nb[0]
                if nb1 is not None:
                    nb2 = list(self.G.neighbors(nb1))
                    if nb2[0] == node: nb2 = nb2[1]
                    elif nb2[1] == node: nb2 = nb2[0]
                    w1 = self.G.edges[nb0, node]['weight']
                    w2 = self.G.edges[node, nb1]['weight']
                    w3 = self.G.edges[nb1, nb2]['weight']
                    if w1 < self.min_len and w2 < self.min_len and w3 < self.min_len:
                        for nb3 in list(self.G.neighbors(nb0)):
                            w0 = w1 + w2 + w3
                            if nb3 is not node:
                                w0 += self.G.edges[nb0, nb3]['weight']
                                b = self.G.edges[nb0, nb3]['pll_branch']
                                self.G.add_edge(nb2, nb3, weight = w0, pll_branch = b)
                        self.G.remove_nodes_from([node, nb0, nb1])

    def remove_short_cycles(self: Self) -> None:
        # Merge junctions across spurious cycles
        cycles = list(nx.minimum_cycle_basis(self.G))
        for i in range(len(cycles)):
            cycle = cycles[i]
            cycle_n = np.asarray([self.G.has_node(node) for node in cycle]).sum()
            if len(cycle) == cycle_n:
                w_cycle = self.G.subgraph(cycle).size('weight')
                if w_cycle < self.min_cycle:
                    cycle_node = f'cycle_{i}_{w_cycle}'
                    nb_cycle, nb_cycle_b = [], []
                    for node in cycle:
                        for nb in list(self.G.neighbors(node)):
                            if nb not in cycle:
                                nb_cycle.append(nb)
                                b = self.G.edges[node, nb]['pll_branch']
                                nb_cycle_b.append(b)
                    for i in range(len(nb_cycle)):
                        nb = nb_cycle[i]
                        b = nb_cycle_b[i]
                        self.G.add_edge(nb, cycle_node, weight = (w_cycle / 2), pll_branch = b)
                    self.G.remove_nodes_from(cycle)

    def simplify_and_prune(self: Self) -> None:
        self.simplify_graph_paths()
        self.prune_overhangs()
        self.simplify_graph_paths()

    def simplify_graph_paths(self: Self) -> None:
        # Leave two nodes on paths between junctions
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 2:
                nb = list(self.G.neighbors(node))
                nb_deg = list(self.G.degree(nb))
                if nb_deg[0][1] == 2 and nb_deg[1][1] == 2:
                    e1 = self.G.edges[node, nb[0]]
                    e2 = self.G.edges[node, nb[1]]
                    w12 = e1['weight'] + e2['weight']
                    b = (e1['pll_branch'] and e2['pll_branch'])
                    self.G.add_edge(*nb, weight = w12, pll_branch = b)
                    self.G.remove_node(node)

    def prune_overhangs(self: Self) -> None:
        # Fix single-node overhangs
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 1:
                nb = list(self.G.neighbors(node))[0]
                if self.G.degree(nb) > 2:
                    self.G.remove_node(node)
        # Fix double-node overhangs
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 1:
                nb = list(self.G.neighbors(node))[0]
                if self.G.degree(nb) == 2:
                    nb1 = list(self.G.neighbors(nb))
                    if nb1[0] == node: nb1 = nb1[1]
                    elif nb1[1] == node: nb1 = nb1[0]
                    if self.G.degree(nb1) > 2:
                        self.G.remove_nodes_from([node, nb])

    def graph_statistics(self: Self) -> None:
        self.cycles = list(nx.minimum_cycle_basis(self.G))
        self.cycles_k = list(map(lambda p: int(len(p) / 3), self.cycles))
        self.cycles_len = self.get_lengths(self.cycles)
        self.cycles_area = self.get_areas(self.cycles)
        self.ends = self.get_free_ends()
        self.ends_len = self.get_lengths(self.ends)
        self.ends_curv = self.get_curvature(self.ends)
        self.ends_z = self.get_depths([x[0] for x in self.ends])
        self.nodes = [node for node, deg in list(self.G.degree()) if deg > 2]
        self.nodes_deg = [deg for _, deg in list(self.G.degree(self.nodes))]
        self.nodes_z = self.get_depths(self.nodes)

    def get_free_ends(self: Self) -> list:
        free_ends = []
        for node in list(self.G.nodes()):
            if self.G.degree(node) == 1:
                nb = list(self.G.neighbors(node))[0]
                nb1 = list(self.G.neighbors(nb))
                if nb1[0] == node: nb1 = nb1[1]
                elif nb1[1] == node: nb1 = nb1[0]
                nb2 = list(self.G.neighbors(nb1))
                if nb2[0] == nb: nb2 = nb2[1]
                elif nb2[1] == nb: nb2 = nb2[0]
                if not self.G.edges[node, nb]['pll_branch']:
                    free_ends.append([node, nb, nb1, nb2])
        return free_ends
    
    def get_lengths(self: Self,
                    nodes: list
                    ) -> list:
        return list(map(lambda H: self.G.subgraph(H).size('weight'), nodes))

    def get_areas(self: Self,
                  cycles: list,
                  min_area: float = 2.
                  ) -> list:
        areas = []
        scaling = np.array([self.z_res, self.xy_res, self.xy_res])
        for cycle in cycles:
            H = self.G.subgraph(cycle)
            nodes = [cycle[0]]
            nb = list(H.neighbors(nodes[0]))[0]
            nodes.append(nb)
            while nodes[-1] != nodes[0]:
                nb = list(H.neighbors(nodes[-1]))
                if nb[0] == nodes[-2]: nb = nb[1]
                elif nb[1] == nodes[-2]: nb = nb[0]
                nodes.append(nb)
            nodes = nodes[:-1]
            try:
                nodes = list(map(lambda n: [float(x) for x in n[1:-1].split(' ')], nodes))
                nodes = np.array(nodes) * scaling[None, :]
                nodes_area = np.zeros((3,))
                for k in range(1, nodes.shape[0] - 1):
                    n0_nk = (nodes[k] - nodes[0])
                    n0_nk1 = (nodes[k + 1] - nodes[0])
                    nodes_area += (.5 * np.cross(n0_nk, n0_nk1))
                nodes_area = np.linalg.norm(nodes_area)
                if nodes_area >= min_area:
                    areas.append(nodes_area)
                else:
                    areas.append(np.nan)
            except:
                areas.append(np.nan)
        return areas
    
    def get_curvature(self: Self,
                      ends: list,
                      sm: float = .3,
                      n_spl: int = 1001
                      ) -> list:
        curvature = []
        scaling = np.array([self.z_res, self.xy_res, self.xy_res])
        for nodes in ends:
            nodes_path = [nodes[0]]
            nb = list(self.H.neighbors(nodes[0]))[0]
            full_path = False
            while not full_path:
                nodes_path.append(nb)
                if self.H.degree(nb) == 2:
                    nb_new = list(self.H.neighbors(nb))
                    if nb_new[0] == nodes_path[-2]: nb = nb_new[1]
                    elif nb_new[1] == nodes_path[-2]: nb = nb_new[0]
                else:
                    if nb == nodes[-1]: full_path = True
                    break
            if full_path:
                try:
                    nodes_path = np.array(list(map(lambda n: [float(x) for x in n[1:-1].split(' ')], nodes_path))) * scaling[None, :]
                    arclen_cum = np.concatenate(([0], np.sqrt((np.diff(nodes_path, axis = 0) ** 2).sum(axis = 1)).cumsum()))
                    tck, _ = splprep([nodes_path[:,0], nodes_path[:,1], nodes_path[:,2]], u = arclen_cum, s = sm)
                    ss = np.linspace(0, arclen_cum.max(), n_spl)
                    zp, yp, xp = splev(ss, tck, der = 1)
                    zpp, ypp, xpp = splev(ss, tck, der = 2)
                    zppp, yppp, xppp = splev(ss, tck, der = 3)
                    curv = ((xp * ypp * zppp) + (yp * zpp * xppp) + (zp * xpp * yppp) - (xp * zpp * yppp) - (yp * xpp * zppp) - (zp * ypp * xppp)) / ((xp ** 2 + yp ** 2 + zp ** 2) ** (3/2))
                    abs_tot_curv = abs(simpson(curv, ss))
                    curvature.append(abs_tot_curv)
                except:
                    curvature.append(np.nan)
            else:
                curvature.append(np.nan)
        return curvature

    def get_depths(self: Self,
                   nodes: list
                   ) -> list:
        scaling = np.array([self.z_res, self.xy_res, self.xy_res])
        pll_branch_node, zyx0 = None, None
        for e in list(self.G.edges()):
            if self.G.edges[e]['pll_branch']:
                if self.G.degree(e[0]) > 2: pll_branch_node = e[0]
                elif self.G.degree(e[1]) > 2: pll_branch_node = e[1]
            if pll_branch_node is not None:
                try: zyx0 = np.array([float(x) for x in pll_branch_node[1:-1].split(' ')])
                except: pass
                break
        depths = [None] * len(nodes)
        for i in range(len(nodes)):
            try:
                node_zyx = np.array([float(x) for x in nodes[i][1:-1].split(' ')])
                u = (zyx0 - node_zyx) * scaling
                v = (zyx0 - self.apex_zyx) * scaling
                depths[i] = (np.dot(u, v) / np.linalg.norm(v)) / np.linalg.norm(v)
            except: 
                depths[i] = np.nan
        return depths

def aggregate_results(res_dir: str) -> Tuple[Union[pd.DataFrame, Tuple[list]]]:
    fn = sorted(list(map(str, pathlib.Path(f'SNT_Excel/{res_dir}').glob('*/*/*.TIF.csv'))))
    df = pd.DataFrame(index = fn, columns = ['n_cycles', 'n_ends', 'n_nodes'])
    cycles_k, cycles_len, cycles_area, cycles_fn = [], [], [], []
    ends_len, ends_curv, ends_z, ends_fn = [], [], [], []
    nodes_deg, nodes_z, nodes_fn = [], [], []
    for ix in df.index:
        print(ix)
        df_ix = pd.read_csv(ix)
        exp = SNT_graph(df_ix)
        df.loc[ix] = [len(exp.cycles), len(exp.ends), len(exp.nodes)]
        cycles_k.extend(exp.cycles_k)
        cycles_len.extend(exp.cycles_len)
        cycles_area.extend(exp.cycles_area)
        cycles_fn.extend([ix] * len(exp.cycles))
        ends_len.extend(exp.ends_len)
        ends_curv.extend(exp.ends_curv)
        ends_z.extend(exp.ends_z)
        ends_fn.extend([ix] * len(exp.ends))
        nodes_deg.extend(exp.nodes_deg)
        nodes_z.extend(exp.nodes_z)
        nodes_fn.extend([ix] * len(exp.nodes))
    return df, (cycles_k, cycles_len, cycles_area, cycles_fn, 
                ends_len, ends_curv, ends_z, ends_fn, 
                nodes_deg, nodes_z, nodes_fn)

_2dpf_wt  = aggregate_results('wt/2dpf')
_3dpf_wt = aggregate_results('wt/3dpf')
_4dpf_wt  = aggregate_results('wt/4dpf')
_2dpf_mut  = aggregate_results('mut/2dpf')
_3dpf_mut  = aggregate_results('mut/3dpf')
_4dpf_mut  = aggregate_results('mut/4dpf')
hc_num = pd.read_csv('SNT_Excel/hc_numbers.csv')

df = pd.concat((_2dpf_wt[0], _3dpf_wt[0], _4dpf_wt[0], _2dpf_mut[0], _3dpf_mut[0], _4dpf_mut[0]), 
               keys = ['2dpf_wt', '3dpf_wt', '4dpf_wt', '2dpf_mut', '3dpf_mut', '4dpf_mut']
               ).reset_index(level = 0)
df.rename({'level_0' : 'condition'}, axis = 'columns', inplace = True)
df['hc_number'] = hc_num['HC Number'].values

#%%
# Node Degree vs. Relative Height
ylim = np.array([-.1, 1.])
wt_nodes_deg = np.array(_4dpf_wt[1][8])
wt_nodes_z = np.array(_4dpf_wt[1][9])
wt_hc_num = df['hc_number'].loc[_4dpf_wt[1][10]].values
n_wt = np.where(df['condition'] == '4dpf_wt')[0].size
ix = np.where(~np.isnan(wt_nodes_z) & (wt_nodes_z >= ylim[0]) & (wt_nodes_z <= ylim[1]))[0]
wt_nodes_deg = wt_nodes_deg[ix]
wt_nodes_z = wt_nodes_z[ix]
wt_hc_num = wt_hc_num[ix]
yy = np.arange(ylim[0], ylim[1] + .05, .1)
wt_nodes_z_binned = np.digitize(wt_nodes_z, yy) - 1
xlim = np.array([3, 7])
xx = np.arange(xlim[0], xlim[1] + 1)
wt = np.zeros((yy.size - 1, xx.size))
for i in range(wt_nodes_deg.size):
    wt[wt_nodes_z_binned[i], wt_nodes_deg[i] - xlim[0]] += (10 / wt_hc_num[i])
wt = np.log10(50 * (wt / n_wt) + 1)   # scale for plotting

mut_nodes_deg = np.array(_4dpf_mut[1][8])
mut_nodes_z = np.array(_4dpf_mut[1][9])
mut_hc_num = df['hc_number'].loc[_4dpf_mut[1][10]].values
n_mut = np.where(df['condition'] == '4dpf_mut')[0].size
ix = np.where(~np.isnan(mut_nodes_z) & (mut_nodes_z >= ylim[0]) & (mut_nodes_z <= ylim[1]))[0]
mut_nodes_deg = mut_nodes_deg[ix]
mut_nodes_z = mut_nodes_z[ix]
mut_hc_num = mut_hc_num[ix]
mut_nodes_z_binned = np.digitize(mut_nodes_z, yy) - 1
mut = np.zeros((yy.size - 1, xx.size))
for i in range(mut_nodes_deg.size):
    mut[mut_nodes_z_binned[i], mut_nodes_deg[i] - xlim[0]] += (10 / mut_hc_num[i])
mut = np.log10(50 * (mut / n_mut) + 1)

fs1, fs2, fs3 = 12, 11, 10
vmax = max([wt.max(), mut.max()])
jet = cm.get_cmap('jet', 256)
colors = jet(np.linspace(0, 1, 256))
cmap = ListedColormap(colors[30 : -30])

fig, ax = plt.subplots(1, 2, figsize = (3.5, 3.5))
hm = ax[0].imshow(wt, origin = 'lower', cmap = cmap, vmax = vmax)
ax[0].set_title('WT', size = fs1)
ax[0].set_xticks(xx - xlim[0], labels = xx, fontsize = fs3)
yticks = (np.arange(yy.size) - .5)
ax[0].set_yticks(yticks, labels = (yy * 100).astype(int), fontsize = fs3)
ax[0].set_ylabel('Relative Height (%)', fontsize = fs2)

ax[1].imshow(mut, origin = 'lower', cmap = cmap, vmax = vmax)
ax[1].set_title('Mutant', size = fs1)
ax[1].set_xticks(xx - xlim[0], labels = xx, fontsize = fs3)
ax[1].set_yticks([])
fig.text(.5, 0, 'Node Degree', ha = 'center', va = 'center', fontsize = fs2)

cax = fig.add_axes([.925, .12, 0.05, .755])
cb = fig.colorbar(hm, cax = cax, ticks = [0, vmax])
cb.ax.set_yticklabels(['0', 'Max'], fontsize = fs3)
fig.text(1.015, .5, 'Relative Number of Nodes', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/degree_vs_height_4dpf.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

pvals = np.zeros(xx.size - 1)
for i in range(xx.size - 1):
    z_wt = wt_nodes_z[wt_nodes_deg == xx[i]]
    z_mut = mut_nodes_z[mut_nodes_deg == xx[i]]
    pvals[i] = kstest(z_wt, z_mut, alternative = 'less').pvalue
print(pvals)

#%%
# Nodes Per Cycle
ylim = np.array([1, 7])
wt_pmf = np.bincount(_4dpf_wt[1][0], minlength = ylim[1])[ylim[0] : (ylim[1] + 1)]
wt_pmf = (wt_pmf / wt_pmf.sum()).reshape(-1, 1) * 100
mut_pmf = np.bincount(_4dpf_mut[1][0], minlength = ylim[1])[ylim[0] : (ylim[1] + 1)]
mut_pmf = (mut_pmf / mut_pmf.sum()).reshape(-1, 1) * 100

fig, ax = plt.subplots(1, 2, figsize = (1., 3.5))
vmin, vmax = 0, max([wt_pmf.max(), mut_pmf.max()]) #36
hm = ax[0].imshow(wt_pmf, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].set_title('WT', size = fs2)
ax[0].set_xticks([])
yy = np.arange(ylim[0], ylim[1] + 1)
ax[0].set_yticks(yy - ylim[0], labels = yy, size = fs3)
ax[0].set_ylabel('Nodes Per Cycle', fontsize = fs2)

ax[1].imshow(mut_pmf, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_title('Mutant', size = fs2)
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_ylim(-.5, ylim[1] - .5)

cax = fig.add_axes([1.05, .15, .175, .685])
cb = fig.colorbar(hm, cax = cax)
fig.text(1.575, .5, '% Total Cycles', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/nodes_per_cycle_4dpf.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

print(kstest(_4dpf_wt[1][0], _4dpf_mut[1][0], alternative = 'less').pvalue)

#%%
# Relative Number of Cycles vs. Free Ends
upper, fit_upper = 12, 8
xx = np.arange(0, upper + 1)
yy = np.arange(0, upper + 1)
wt = np.zeros((yy.size, xx.size))
mut = np.zeros((yy.size, xx.size))
n_wt = np.where(df['condition'] == '4dpf_wt')[0].size
n_mut = np.where(df['condition'] == '4dpf_mut')[0].size
ii_fit, jj_fit = [], []
for ix in df.index:
    ii, jj = df[['n_ends', 'n_cycles']].loc[ix].values
    ii = 10 * (ii / df['hc_number'].loc[ix])   # scaled
    jj = 10 * (jj / df['hc_number'].loc[ix])
    ii = np.digitize(ii, yy)
    jj = np.digitize(jj, xx)
    if ii <= upper and jj <= upper:
        if df['condition'].loc[ix] == '4dpf_wt':
            wt[ii, jj] += 1
            ii_fit.append(ii)
            jj_fit.append(jj)
        elif df['condition'].loc[ix] == '4dpf_mut':
            mut[ii, jj] += 1
wt /= n_wt
mut /= n_mut
ii_fit = np.array(ii_fit)
jj_fit = np.array(jj_fit)
ix_fit = (ii_fit < fit_upper) & (jj_fit < fit_upper)
ii_fit = ii_fit[ix_fit]
jj_fit = jj_fit[ix_fit]
jj_fit = jj_fit.reshape(-1, 1)
reg = LinearRegression().fit(jj_fit, ii_fit)
_jj_fit = np.arange(1, 11).reshape(-1, 1)
_ii_fit = reg.predict(_jj_fit)
r2, pval = pearsonr(jj_fit.flatten(), ii_fit)

fig, ax = plt.subplots(1, 2, figsize = (6., 3.))
vmin, vmax = 0, max([wt.max(), mut.max()])
hm = ax[0].imshow(wt, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].plot(_jj_fit, _ii_fit, 'r--', label = r'$R^2 = $' + f'{r2.round(2)}')
ax[0].set_title('WT', size = fs1)
ticks = np.linspace(0, upper, 6)
ticklabels = np.arange(0, 1.1, .2).round(1)
ax[0].set_xticks(ticks, labels = ticklabels, fontsize = fs3)
ax[0].set_yticks(ticks, labels = ticklabels, fontsize = fs3)
ax[0].set_ylabel('Rel. Number of Free Ends', fontsize = fs2)
ax[0].legend(frameon = False, fontsize = fs3)

ax[1].imshow(mut, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_title('Mutant', size = fs1)
ax[1].set_xticks(ticks, labels = ticklabels, fontsize = fs3)
ax[1].set_yticks([])
fig.text(.5, 0, 'Relative Number of Cycles', ha = 'center', va = 'center', fontsize = fs2)

cax = fig.add_axes([.925, .15, .025, .7])
cb = fig.colorbar(hm, cax = cax, ticks = [vmin, vmax])
cb.ax.set_yticklabels(['0', 'Max'], fontsize = fs3)
fig.text(.975, .5, 'Normalized Density', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/num_cycles_vs_ends_4dpf.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

print(pval)

#%%
# Contour Length: % Free Ends vs. Cycles
ylim = np.array([0, 70])
cycles_len = np.array(_4dpf_wt[1][1])
ends_len = np.array(_4dpf_wt[1][4])
ix_cycles = np.where(cycles_len < ylim[1])[0]
cycles_len = cycles_len[ix_cycles]
ix_ends = np.where(ends_len < ylim[1])[0]
ends_len = ends_len[ix_ends]
yy = np.arange(ylim[0], ylim[1] + 1, 10)
cycles_len_binned = np.digitize(cycles_len, yy) - 1
ends_len_binned = np.digitize(ends_len, yy) - 1
len_arr = np.zeros((yy.size - 1, 2))
for i in range(ends_len.size):
    len_arr[ends_len_binned[i], 0] += 1
for i in range(cycles_len.size):
    len_arr[cycles_len_binned[i], 1] += 1
len_arr = (len_arr / len_arr.sum(axis = 1)[:,None]) * 100
len_arr[np.isnan(len_arr)] = 0
ends_arr = len_arr[:, 0].reshape(-1, 1)
cycles_arr = len_arr[:, 1].reshape(-1, 1)

fig, ax = plt.subplots(1, 2, figsize = (1., 3.5))
vmin, vmax = 0, len_arr.max()
hm = ax[0].imshow(ends_arr, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].set_xticks([0], labels = ['Ends'], size = fs3)
yticks = np.arange(ylim[0], (ylim[1] / 10) + 1) - .5
ax[0].set_yticks(yticks, labels = yy, size = fs3)
ax[0].set_ylabel(r'Contour Length ($\mu m$)', fontsize = fs2)
fig.text(.5, .89, 'WT', ha = 'center', va = 'center', size = fs1)

ax[1].imshow(cycles_arr, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_xticks([0], labels = ['Cycles'], size = fs3)
ax[1].set_yticks([])

cax = fig.add_axes([1., .155, .175, .68])
cb = fig.colorbar(hm, cax = cax)
fig.text(1.5, .5, '% Total Paths', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/len_cycles_vs_ends_4dpf_wt.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

# wt_freq = len_arr[-3:-1, :].sum(0)
# wt_freq = (wt_freq / wt_freq.sum()) * 100
# mut_freq = len_arr[-3:-1, :].sum(0)
# mut_freq = (mut_freq / mut_freq.sum()) * 100
# print(chisquare(mut_freq, wt_freq).pvalue)

#%%
# Relative Height of Bare Terminals
ylim = np.array([0, 1])
wt_ends_z = np.array(_3dpf_wt[1][6])
ix = np.where(~np.isnan(wt_ends_z) & (wt_ends_z >= ylim[0]) & (wt_ends_z <= ylim[1]))[0]
wt_ends_z = wt_ends_z[ix]
yy = np.arange(ylim[0], ylim[1] + .1, .2)
wt_ends_z_binned = np.digitize(wt_ends_z, yy) - 1
wt = np.zeros((yy.size - 1, 1))
for i in range(wt_ends_z.size):
    wt[wt_ends_z_binned[i]] += 1
wt = (wt / wt.sum()) * 100

mut_ends_z = np.array(_3dpf_mut[1][6])
ix = np.where(~np.isnan(mut_ends_z) & (mut_ends_z >= ylim[0]) & (mut_ends_z <= ylim[1]))[0]
mut_ends_z = mut_ends_z[ix]
mut_ends_z_binned = np.digitize(mut_ends_z, yy) - 1
mut = np.zeros((yy.size - 1, 1))
for i in range(mut_ends_z.size):
    mut[mut_ends_z_binned[i]] += 1
mut = (mut / mut.sum()) * 100

fig, ax = plt.subplots(1, 2, figsize = (.95, 3.75))
vmin, vmax = 0, 31
hm = ax[0].imshow(wt, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].set_title('WT', size = fs3)
ax[0].set_xticks([])
yticks = np.arange(0, yy.size) - .5
ax[0].set_yticks(yticks, labels = (yy *  100).astype(int), size = fs3)
ax[0].set_ylabel('Relative Height (%)', fontsize = fs3)

ax[1].imshow(mut, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_title('Mutant', size = fs3)
ax[1].set_xticks([])
ax[1].set_yticks([])

cax = fig.add_axes([1.05, .275, .175, .44])
cb = fig.colorbar(hm, cax = cax)
fig.text(1.575, .5, '% Terminals', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/terminals_height_3dpf.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

print(kstest(wt_ends_z, mut_ends_z, alternative = 'less').pvalue)

#%%
# Free Ends: Curvature vs. Length
ylim = np.array([0, 30])
curv = np.array(_4dpf_wt[1][5])
ends_len = np.array(_4dpf_wt[1][4])
curv_thresh = (np.pi / 6)
idx_low = (curv < curv_thresh) & (ends_len <= ylim[1])
idx_high = (curv > curv_thresh) & (ends_len <= ylim[1])
len_low = ends_len[idx_low]
len_high = ends_len[idx_high]
yy = np.arange(ylim[0], ylim[1] + 1, 5)
len_low_binned = np.digitize(len_low, yy) - 1
len_high_binned = np.digitize(len_high, yy) - 1
len_arr = np.zeros((yy.size - 1, 2))
for i in range(len_low.size):
    len_arr[len_low_binned[i], 0] += 1
for i in range(len_high.size):
    len_arr[len_high_binned[i], 1] += 1
len_arr = (len_arr / len_arr.sum(axis = None)) * 100
len_arr[np.isnan(len_arr)] = 0
len_low_arr = len_arr[:, 0].reshape(-1, 1)
len_high_arr = len_arr[:, 1].reshape(-1, 1)

fig, ax = plt.subplots(1, 2, figsize = (1., 3.5))
vmin, vmax = 0, 29 #33
hm = ax[0].imshow(len_low_arr, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].set_xticks([0], labels = ['Low'], size = fs3)
yticks = np.arange(ylim[0], (ylim[1] / 5) + 1) - .5
ax[0].set_yticks(yticks, labels = yy, size = fs3)
ax[0].set_ylabel(r'Contour Length ($\mu m$)', fontsize = fs2)
fig.text(.5, .085, 'Curvature', ha = 'center', va = 'center', size = fs3)
fig.text(.5, .85, 'WT', ha = 'center', va = 'center', size = fs1)

ax[1].imshow(len_high_arr, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_xticks([0], labels = ['High'], size = fs3)
ax[1].set_yticks([])

cax = fig.add_axes([.98, .2, .16, .59])
cb = fig.colorbar(hm, cax = cax)
fig.text(1.5, .5, '% Occurences', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/curvature_vs_len_4dpf_wt.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

# wt_len = len_high
# mut_len = len_high
# print(kstest(wt_len, mut_len, alternative = 'less').pvalue)

#%%
# Contour Lengths of Loops
ylim = np.array([0, 70])
cycles_len_wt = np.array(_4dpf_wt[1][1])
cycles_len_mut = np.array(_4dpf_mut[1][1])
ix_wt = np.where(cycles_len_wt < ylim[1])[0]
cycles_len_wt = cycles_len_wt[ix_wt]
ix_mut = np.where(cycles_len_mut < ylim[1])[0]
cycles_len_mut = cycles_len_mut[ix_mut]
yy = np.arange(ylim[0], ylim[1] + 1, 10)
len_wt_binned = np.digitize(cycles_len_wt, yy) - 1
len_mut_binned = np.digitize(cycles_len_mut, yy) - 1
len_arr = np.zeros((yy.size - 1, 2))
for i in range(cycles_len_wt.size):
    len_arr[len_wt_binned[i], 0] += 1
for i in range(cycles_len_mut.size):
    len_arr[len_mut_binned[i], 1] += 1
len_arr = (len_arr / len_arr.sum(axis = 0)[None,:]) * 100
len_arr[np.isnan(len_arr)] = 0
len_arr_wt = len_arr[:, 0].reshape(-1, 1)
len_arr_mut = len_arr[:, 1].reshape(-1, 1)

fig, ax = plt.subplots(1, 2, figsize = (1., 3.5))
vmin, vmax = 0, 45
hm = ax[0].imshow(len_arr_wt, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[0].set_title('WT', size = fs2)
ax[0].set_xticks([])
yticks = np.arange(ylim[0], (ylim[1] / 10) + 1) - .5
ax[0].set_yticks(yticks, labels = yy, size = fs3)
ax[0].set_ylabel(r'Contour Length ($\mu m$)', fontsize = fs2)

ax[1].imshow(len_arr_mut, origin = 'lower', cmap = 'Blues', vmin = vmin, vmax = vmax)
ax[1].set_title('Mutant', size = fs2)
ax[1].set_xticks([])
ax[1].set_yticks([])

cax = fig.add_axes([1.05, .15, .175, .685])
cb = fig.colorbar(hm, cax = cax)
fig.text(1.575, .5, '% Cycles', rotation = -90, va = 'center', fontsize = fs2)
# plt.savefig('snt_graph_figures/len_cycles_4dpf.pdf', bbox_inches = 'tight', dpi = 600)
plt.show()

print(kstest(_4dpf_wt[1][1], _4dpf_mut[1][1], alternative = 'less').pvalue)