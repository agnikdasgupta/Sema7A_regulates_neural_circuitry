import os
import sys
from tkinter import Tk, filedialog
import matplotlib.pyplot as plt
import numpy as np
import pandas
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
from natsort import natsorted

import itertools
import matplotlib.cm as cm

import seaborn as sns


folder_path = 'E:/220625-Sema71-Arborization/WT/2 dpf complete HQ/Pooled top_and_bot - 17'

csv_name = (natsorted(j for j in os.listdir(f'{folder_path}') if (f'-top_and_bot' in j) and (f'.csv' in j)))

x = []
y = []
z = []
x_offset = []
y_offset = []
z_offset = []

df = []
for i in range(len(csv_name)):
    df.append(pd.read_csv(f'{folder_path}/{csv_name[i]}'))
    x.append(df[i]["X"].values)
    y.append(df[i]["Y"].values)
    z.append(df[i]["Z"].values)

    x_offset.append(x[i][0])
    y_offset.append(y[i][0])
    z_offset.append(z[i][-1])
    
    x[i][0] = x[i][-2]
    y[i][0] = y[i][-2]
    z[i][0] = z[i][-2]
 

    z[i] = [j-z_offset[i] for j in z[i]]
    z[i] = [-j for j in z[i]]
    z[i] = [round(0.2*j,4) for j in z[i]]
    x[i] = [j-x_offset[i] for j in x[i]]
    x[i] = [round(0.065*j,4) for j in x[i]]
    y[i] = [j-y_offset[i] for j in y[i]]
    y[i] = [round(0.065*j,4) for j in y[i]]
    
distEOB = []
distO = []
for i in range(len(x)):
    distEOB.append([0]*len(x[i]))
    distO.append([0]*len(x[i]))
    for j in range(len(x[i])):
        distEOB[i][j] = np.sqrt((x[i][j]-x[i][-1])**2+(y[i][j]-y[i][-1])**2+(z[i][j]-z[i][-1])**2)
        distO[i][j] = np.sqrt((x[i][j])**2+(y[i][j])**2+(z[i][j])**2)

distEOBlin = []
for value in distEOB:
    for subvalue in value:
        distEOBlin.append(subvalue)
distOlin = []
for value in distO:
    for subvalue in value:
        distOlin.append(subvalue)

x[i][-1] = x[i][-2]
y[i][-1] = y[i][-2]
z[i][-1] = z[i][-2]

# fig = plt.figure(figsize=(9.6,9.6))
fig = plt.figure(figsize=(11.9,9.6))
ax = fig.add_subplot(111, projection='3d')

ax.set_xticks([-30,-20,-10,0,10,20,30])
ax.set_yticks([-30,-20,-10,0,10,20,30])
ax.set_zticks([0,10,20,30,40,50,60])
ax.set_zticks([0,10,20,30,40])

# ax.set_facecolor("0.50")

ax.axes.set_xlim3d(left=-30, right=30) 
ax.axes.set_ylim3d(bottom=-30, top=30) 
# ax.axes.set_xlim3d(left=-10, right=10)
# ax.axes.set_ylim3d(bottom=-10, top=10)
ax.axes.set_zlim3d(bottom=0, top=45)

pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='jet',s=1,vmin=0,vmax=25)
# pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='viridis',s=1,vmin=0,vmax=1)
cbar = plt.colorbar(pnt3d)
cbar.set_label("Distance from PLL branch (µm)")
for i in range(len(csv_name)):
# for i in [-3]:
    # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1,vmin=0,vmax=25)
    ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1)
    # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='viridis',s=1,vmin=0,vmax=1)
    

# ax.view_init(90, 90)
# ax.view_init(5, 90)
ax.view_init(15, 90)

# for i in range(360):
#     ax.azim += 1
#     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-17-2dpf-movie_f-WT-{i}.png')

plt.ion()
# plt.savefig("E:/Workstation/Skeletonization_Excel/plot1.pdf")

plt.show()






# MUTANT

folder_path2 = 'E:/220625-Sema71-Arborization/sa24691/2 dpf complete HQ/Pooled top_and_bot - 17'

csv_name2 = (natsorted(j for j in os.listdir(f'{folder_path2}') if (f'-top_and_bot' in j) and (f'.csv' in j)))
x2 = []
y2 = []
z2 = []
x_offset2 = []
y_offset2 = []
z_offset2 = []

df2 = []
for i in range(len(csv_name2)):
    df2.append(pd.read_csv(f'{folder_path2}/{csv_name2[i]}'))
    x2.append(df2[i]["X"].values)
    y2.append(df2[i]["Y"].values)
    z2.append(df2[i]["Z"].values)

    x_offset2.append(x2[i][0])
    y_offset2.append(y2[i][0])
    z_offset2.append(z2[i][-1])
    
    x2[i][0] = x2[i][-2]
    y2[i][0] = y2[i][-2]
    z2[i][0] = z2[i][-2]


    z2[i] = [j-z_offset2[i] for j in z2[i]]
    z2[i] = [-j for j in z2[i]]
    z2[i] = [round(0.2*j,4) for j in z2[i]]
    x2[i] = [j-x_offset2[i] for j in x2[i]]
    x2[i] = [round(0.065*j,4) for j in x2[i]]
    y2[i] = [j-y_offset2[i] for j in y2[i]]
    y2[i] = [round(0.065*j,4) for j in y2[i]]
    
distEOB2 = []
distO2 = []
for i in range(len(x2)):
    distEOB2.append([0]*len(x2[i]))
    distO2.append([0]*len(x2[i]))
    for j in range(len(x2[i])):
        distEOB2[i][j] = np.sqrt((x2[i][j]-x2[i][-1])**2+(y2[i][j]-y2[i][-1])**2+(z2[i][j]-z2[i][-1])**2)
        distO2[i][j] = np.sqrt((x2[i][j])**2+(y2[i][j])**2+(z2[i][j])**2)

distEOBlin2 = []
for value in distEOB2:
    for subvalue in value:
        distEOBlin2.append(subvalue)
distOlin2 = []
for value in distO2:
    for subvalue in value:
        distOlin2.append(subvalue)

x2[i][-1] = x2[i][-2]
y2[i][-1] = y2[i][-2]
z2[i][-1] = z2[i][-2]


# fig = plt.figure(figsize=(9.6,9.6))
fig = plt.figure(figsize=(11.9,9.6))
ax = fig.add_subplot(111, projection='3d')

ax.set_xticks([-30,-20,-10,0,10,20,30])
ax.set_yticks([-30,-20,-10,0,10,20,30])
ax.set_zticks([0,10,20,30,40,50,60])
ax.set_zticks([0,10,20,30,40])

ax.axes.set_xlim3d(left=-30, right=30) 
ax.axes.set_ylim3d(bottom=-30, top=30)
# ax.axes.set_xlim3d(left=-10, right=10)
# ax.axes.set_ylim3d(bottom=-10, top=10)
ax.axes.set_zlim3d(bottom=0, top=45)

# ax.set_facecolor("0.75")

pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='jet',s=1,vmin=0,vmax=25)
# pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='viridis',s=1,vmin=0,vmax=1)

cbar2 = plt.colorbar(pnt3d2)
cbar2.set_label("Distance from PLL branch (µm)")
for i in range(len(csv_name2)):
# for i in [-3]:
    # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1,vmin=0,vmax=25)
    ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1)
    # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='viridis',s=1,vmin=0,vmax=1)

    
# ax.view_init(90, 90)
# ax.view_init(5, 90)
ax.view_init(15, 90)


# for i in range(360):
#     ax.azim += 1
#     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-17-2dpf-movie_f-MT-{i}.png')

plt.ion()
# plt.savefig("E:/Workstation/Skeletonization_Excel/plot2.pdf")

plt.show()





print(np.average(distEOBlin))
print(np.average(distEOBlin2))
print(np.average(distOlin))
print(np.average(distOlin2))



# folder_path = 'E:/220625-Sema71-Arborization/WT/3 dpf complete HQ/Pooled top_and_bot - 29'

# csv_name = (natsorted(j for j in os.listdir(f'{folder_path}') if (f'-top_and_bot' in j) and (f'.csv' in j)))

# x = []
# y = []
# z = []
# x_offset = []
# y_offset = []
# z_offset = []

# df = []
# for i in range(len(csv_name)):
#     df.append(pd.read_csv(f'{folder_path}/{csv_name[i]}'))
#     x.append(df[i]["X"].values)
#     y.append(df[i]["Y"].values)
#     z.append(df[i]["Z"].values)

#     x_offset.append(x[i][0])
#     y_offset.append(y[i][0])
#     z_offset.append(z[i][-1])
    
#     x[i][0] = x[i][-2]
#     y[i][0] = y[i][-2]
#     z[i][0] = z[i][-2]
 

#     z[i] = [j-z_offset[i] for j in z[i]]
#     z[i] = [-j for j in z[i]]
#     z[i] = [round(0.2*j,4) for j in z[i]]
#     x[i] = [j-x_offset[i] for j in x[i]]
#     x[i] = [round(0.065*j,4) for j in x[i]]
#     y[i] = [j-y_offset[i] for j in y[i]]
#     y[i] = [round(0.065*j,4) for j in y[i]]
    
# distEOB = []
# distO = []
# for i in range(len(x)):
#     distEOB.append([0]*len(x[i]))
#     distO.append([0]*len(x[i]))
#     for j in range(len(x[i])):
#         distEOB[i][j] = np.sqrt((x[i][j]-x[i][-1])**2+(y[i][j]-y[i][-1])**2+(z[i][j]-z[i][-1])**2)
#         distO[i][j] = np.sqrt((x[i][j])**2+(y[i][j])**2+(z[i][j])**2)

# distEOBlin = []
# for value in distEOB:
#     for subvalue in value:
#         distEOBlin.append(subvalue)
# distOlin = []
# for value in distO:
#     for subvalue in value:
#         distOlin.append(subvalue)

# x[i][-1] = x[i][-2]
# y[i][-1] = y[i][-2]
# z[i][-1] = z[i][-2]

# # fig = plt.figure(figsize=(9.6,9.6))
# fig = plt.figure(figsize=(11.9,9.6))
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xticks([-30,-20,-10,0,10,20,30])
# ax.set_yticks([-30,-20,-10,0,10,20,30])
# ax.set_zticks([0,10,20,30,40,50,60])
# ax.set_zticks([0,10,20,30,40])

# # ax.set_facecolor("0.50")

# ax.axes.set_xlim3d(left=-30, right=30) 
# ax.axes.set_ylim3d(bottom=-30, top=30) 
# # ax.axes.set_xlim3d(left=-10, right=10)
# # ax.axes.set_ylim3d(bottom=-10, top=10)
# ax.axes.set_zlim3d(bottom=0, top=45)

# pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='jet',s=1,vmin=0,vmax=25)
# # pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='viridis',s=1,vmin=0,vmax=1)
# cbar = plt.colorbar(pnt3d)
# cbar.set_label("Distance from PLL branch (µm)")
# for i in range(len(csv_name)):
# # for i in [-3]:
#     # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1,vmin=0,vmax=25)
#     ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1)
#     # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='viridis',s=1,vmin=0,vmax=1)
    

# # ax.view_init(90, 90)
# # ax.view_init(5, 90)
# ax.view_init(15, 90)

# # for i in range(360):
# #     ax.azim += 1
# #     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-29-3dpf-movie_f-WT-{i}.png')

# plt.ion()
# # plt.savefig("E:/Workstation/Skeletonization_Excel/plot1.pdf")

# plt.show()






# # MUTANT

# folder_path2 = 'E:/220625-Sema71-Arborization/sa24691/3 dpf complete HQ/Pooled top_and_bot - 29'

# csv_name2 = (natsorted(j for j in os.listdir(f'{folder_path2}') if (f'-top_and_bot' in j) and (f'.csv' in j)))
# x2 = []
# y2 = []
# z2 = []
# x_offset2 = []
# y_offset2 = []
# z_offset2 = []

# df2 = []
# for i in range(len(csv_name2)):
#     df2.append(pd.read_csv(f'{folder_path2}/{csv_name2[i]}'))
#     x2.append(df2[i]["X"].values)
#     y2.append(df2[i]["Y"].values)
#     z2.append(df2[i]["Z"].values)

#     x_offset2.append(x2[i][0])
#     y_offset2.append(y2[i][0])
#     z_offset2.append(z2[i][-1])
    
#     x2[i][0] = x2[i][-2]
#     y2[i][0] = y2[i][-2]
#     z2[i][0] = z2[i][-2]


#     z2[i] = [j-z_offset2[i] for j in z2[i]]
#     z2[i] = [-j for j in z2[i]]
#     z2[i] = [round(0.2*j,4) for j in z2[i]]
#     x2[i] = [j-x_offset2[i] for j in x2[i]]
#     x2[i] = [round(0.065*j,4) for j in x2[i]]
#     y2[i] = [j-y_offset2[i] for j in y2[i]]
#     y2[i] = [round(0.065*j,4) for j in y2[i]]
    
# distEOB2 = []
# distO2 = []
# for i in range(len(x2)):
#     distEOB2.append([0]*len(x2[i]))
#     distO2.append([0]*len(x2[i]))
#     for j in range(len(x2[i])):
#         distEOB2[i][j] = np.sqrt((x2[i][j]-x2[i][-1])**2+(y2[i][j]-y2[i][-1])**2+(z2[i][j]-z2[i][-1])**2)
#         distO2[i][j] = np.sqrt((x2[i][j])**2+(y2[i][j])**2+(z2[i][j])**2)

# distEOBlin2 = []
# for value in distEOB2:
#     for subvalue in value:
#         distEOBlin2.append(subvalue)
# distOlin2 = []
# for value in distO2:
#     for subvalue in value:
#         distOlin2.append(subvalue)

# x2[i][-1] = x2[i][-2]
# y2[i][-1] = y2[i][-2]
# z2[i][-1] = z2[i][-2]


# # fig = plt.figure(figsize=(9.6,9.6))
# fig = plt.figure(figsize=(11.9,9.6))
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xticks([-30,-20,-10,0,10,20,30])
# ax.set_yticks([-30,-20,-10,0,10,20,30])
# ax.set_zticks([0,10,20,30,40,50,60])
# ax.set_zticks([0,10,20,30,40])

# ax.axes.set_xlim3d(left=-30, right=30) 
# ax.axes.set_ylim3d(bottom=-30, top=30)
# # ax.axes.set_xlim3d(left=-10, right=10)
# # ax.axes.set_ylim3d(bottom=-10, top=10)
# ax.axes.set_zlim3d(bottom=0, top=45)

# # ax.set_facecolor("0.75")

# pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='jet',s=1,vmin=0,vmax=25)
# # pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='viridis',s=1,vmin=0,vmax=1)

# cbar2 = plt.colorbar(pnt3d2)
# cbar2.set_label("Distance from PLL branch (µm)")
# for i in range(len(csv_name2)):
# # for i in [-3]:
#     # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1,vmin=0,vmax=25)
#     ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1)
#     # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='viridis',s=1,vmin=0,vmax=1)

    
# # ax.view_init(90, 90)
# # ax.view_init(5, 90)
# ax.view_init(15, 90)


# # for i in range(360):
# #     ax.azim += 1
# #     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-29-3dpf-movie_f-MT-{i}.png')

# plt.ion()
# # plt.savefig("E:/Workstation/Skeletonization_Excel/plot2.pdf")

# plt.show()








# folder_path = 'E:/220625-Sema71-Arborization/WT/4 dpf complete HQ/Pooled top_and_bot - 27'

# csv_name = (natsorted(j for j in os.listdir(f'{folder_path}') if (f'-top_and_bot' in j) and (f'.csv' in j)))

# x = []
# y = []
# z = []
# x_offset = []
# y_offset = []
# z_offset = []

# df = []
# for i in range(len(csv_name)):
#     df.append(pd.read_csv(f'{folder_path}/{csv_name[i]}'))
#     x.append(df[i]["X"].values)
#     y.append(df[i]["Y"].values)
#     z.append(df[i]["Z"].values)

#     x_offset.append(x[i][0])
#     y_offset.append(y[i][0])
#     z_offset.append(z[i][-1])
    
#     x[i][0] = x[i][-2]
#     y[i][0] = y[i][-2]
#     z[i][0] = z[i][-2]
 

#     z[i] = [j-z_offset[i] for j in z[i]]
#     z[i] = [-j for j in z[i]]
#     z[i] = [round(0.2*j,4) for j in z[i]]
#     x[i] = [j-x_offset[i] for j in x[i]]
#     x[i] = [round(0.065*j,4) for j in x[i]]
#     y[i] = [j-y_offset[i] for j in y[i]]
#     y[i] = [round(0.065*j,4) for j in y[i]]
    
# distEOB = []
# distO = []
# for i in range(len(x)):
#     distEOB.append([0]*len(x[i]))
#     distO.append([0]*len(x[i]))
#     for j in range(len(x[i])):
#         distEOB[i][j] = np.sqrt((x[i][j]-x[i][-1])**2+(y[i][j]-y[i][-1])**2+(z[i][j]-z[i][-1])**2)
#         distO[i][j] = np.sqrt((x[i][j])**2+(y[i][j])**2+(z[i][j])**2)

# distEOBlin = []
# for value in distEOB:
#     for subvalue in value:
#         distEOBlin.append(subvalue)
# distOlin = []
# for value in distO:
#     for subvalue in value:
#         distOlin.append(subvalue)

# x[i][-1] = x[i][-2]
# y[i][-1] = y[i][-2]
# z[i][-1] = z[i][-2]

# # fig = plt.figure(figsize=(9.6,9.6))
# fig = plt.figure(figsize=(11.9,9.6))
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xticks([-30,-20,-10,0,10,20,30])
# ax.set_yticks([-30,-20,-10,0,10,20,30])
# ax.set_zticks([0,10,20,30,40,50,60])
# ax.set_zticks([0,10,20,30,40])

# # ax.set_facecolor("0.50")

# ax.axes.set_xlim3d(left=-30, right=30) 
# ax.axes.set_ylim3d(bottom=-30, top=30) 
# # ax.axes.set_xlim3d(left=-10, right=10)
# # ax.axes.set_ylim3d(bottom=-10, top=10)
# ax.axes.set_zlim3d(bottom=0, top=45)

# pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='jet',s=1,vmin=0,vmax=25)
# # pnt3d = ax.scatter(x[0],y[0],z[0],c=distEOB[0],cmap='viridis',s=1,vmin=0,vmax=1)
# cbar = plt.colorbar(pnt3d)
# cbar.set_label("Distance from PLL branch (µm)")
# for i in range(len(csv_name)):
# # for i in [-3]:
#     # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1,vmin=0,vmax=25)
#     ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='jet',s=1)
#     # ax.scatter(x[i],y[i],z[i],c=distEOB[i],cmap='viridis',s=1,vmin=0,vmax=1)
    

# # ax.view_init(90, 90)
# # ax.view_init(5, 90)
# ax.view_init(15, 90)

# # for i in range(360):
# #     ax.azim += 1
# #     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-27-4dpf-movie_f-WT-{i}.png')

# plt.ion()
# # plt.savefig("E:/Workstation/Skeletonization_Excel/plot1.pdf")

# plt.show()






# # MUTANT

# folder_path2 = 'E:/220625-Sema71-Arborization/sa24691/4 dpf complete HQ/Pooled top_and_bot - 27'

# csv_name2 = (natsorted(j for j in os.listdir(f'{folder_path2}') if (f'-top_and_bot' in j) and (f'.csv' in j)))
# x2 = []
# y2 = []
# z2 = []
# x_offset2 = []
# y_offset2 = []
# z_offset2 = []

# df2 = []
# for i in range(len(csv_name2)):
#     df2.append(pd.read_csv(f'{folder_path2}/{csv_name2[i]}'))
#     x2.append(df2[i]["X"].values)
#     y2.append(df2[i]["Y"].values)
#     z2.append(df2[i]["Z"].values)

#     x_offset2.append(x2[i][0])
#     y_offset2.append(y2[i][0])
#     z_offset2.append(z2[i][-1])
    
#     x2[i][0] = x2[i][-2]
#     y2[i][0] = y2[i][-2]
#     z2[i][0] = z2[i][-2]


#     z2[i] = [j-z_offset2[i] for j in z2[i]]
#     z2[i] = [-j for j in z2[i]]
#     z2[i] = [round(0.2*j,4) for j in z2[i]]
#     x2[i] = [j-x_offset2[i] for j in x2[i]]
#     x2[i] = [round(0.065*j,4) for j in x2[i]]
#     y2[i] = [j-y_offset2[i] for j in y2[i]]
#     y2[i] = [round(0.065*j,4) for j in y2[i]]
    
# distEOB2 = []
# distO2 = []
# for i in range(len(x2)):
#     distEOB2.append([0]*len(x2[i]))
#     distO2.append([0]*len(x2[i]))
#     for j in range(len(x2[i])):
#         distEOB2[i][j] = np.sqrt((x2[i][j]-x2[i][-1])**2+(y2[i][j]-y2[i][-1])**2+(z2[i][j]-z2[i][-1])**2)
#         distO2[i][j] = np.sqrt((x2[i][j])**2+(y2[i][j])**2+(z2[i][j])**2)

# distEOBlin2 = []
# for value in distEOB2:
#     for subvalue in value:
#         distEOBlin2.append(subvalue)
# distOlin2 = []
# for value in distO2:
#     for subvalue in value:
#         distOlin2.append(subvalue)

# x2[i][-1] = x2[i][-2]
# y2[i][-1] = y2[i][-2]
# z2[i][-1] = z2[i][-2]


# # fig = plt.figure(figsize=(9.6,9.6))
# fig = plt.figure(figsize=(11.9,9.6))
# ax = fig.add_subplot(111, projection='3d')

# ax.set_xticks([-30,-20,-10,0,10,20,30])
# ax.set_yticks([-30,-20,-10,0,10,20,30])
# ax.set_zticks([0,10,20,30,40,50,60])
# ax.set_zticks([0,10,20,30,40])

# ax.axes.set_xlim3d(left=-30, right=30) 
# ax.axes.set_ylim3d(bottom=-30, top=30)
# # ax.axes.set_xlim3d(left=-10, right=10)
# # ax.axes.set_ylim3d(bottom=-10, top=10)
# ax.axes.set_zlim3d(bottom=0, top=45)

# # ax.set_facecolor("0.75")

# pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='jet',s=1,vmin=0,vmax=25)
# # pnt3d2 = ax.scatter(x2[0],y2[0],z2[0],c=distEOB2[0],cmap='viridis',s=1,vmin=0,vmax=1)

# cbar2 = plt.colorbar(pnt3d2)
# cbar2.set_label("Distance from PLL branch (µm)")
# for i in range(len(csv_name2)):
# # for i in [-3]:
#     # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1,vmin=0,vmax=25)
#     ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='jet',s=1)
#     # ax.scatter(x2[i],y2[i],z2[i],c=distEOB2[i],cmap='viridis',s=1,vmin=0,vmax=1)

    
# # ax.view_init(90, 90)
# # ax.view_init(5, 90)
# ax.view_init(15, 90)


# # for i in range(360):
# #     ax.azim += 1
# #     plt.savefig(f'E:/Workstation/Skeletonization_Excel/AA-27-4dpf-movie_f-MT-{i}.png')

# plt.ion()
# # plt.savefig("E:/Workstation/Skeletonization_Excel/plot2.pdf")

# plt.show()




# # HIST




hist_min = 0
hist_max = 60
hist_bins = np.linspace(hist_min,hist_max,3601)

hist_bins2 = []
for i in range(3601):
    hist_bins2.append(np.sqrt(i))


counts, bin_edges = np.histogram(distOlin, bins=hist_bins2)
counts2, bin_edges = np.histogram(distOlin2, bins=hist_bins2)
output = open(f'{folder_path}/4dpf_WT_counts.txt','a')
for value in counts:
    output.write(f'{value}\n')
output.close()
output = open(f'{folder_path}/4dpf_sa91_counts.txt','a')
for value in counts2:
    output.write(f'{value}\n')
output.close()


plt.figure()
plt.hist(distOlin,hist_bins2)
# plt.figure()
plt.hist(distOlin2,hist_bins2)

# counts,bins,bars = plt.hist(distOlin,bins=hist_bins)
# disthist = []
# for i in range(len(counts)):
#     counts[i] = counts[i]/((1+2*i)*np.pi)
#     for j in range(counts[i]):
#         disthist.append(i)
    
# counts2,bins,bars = plt.hist(distOlin2,bins=hist_bins)
# for i in range(len(counts)):
#     counts2[i] = counts2[i]/((1+2*i)*np.pi)

plt.figure()


sns.set_style("white")
sns.distplot(distOlin, color="dodgerblue", kde=False, label="Control", hist_kws={'alpha':.3}, bins=hist_bins2)
sns.distplot(distOlin2, color="deeppink", kde=False, label="Mutant", hist_kws={'alpha':.3},bins=hist_bins2)
# sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.yscale('log')

# plt.ylim(0,0.11)
plt.legend()

plt.savefig("E:/Workstation/Skeletonization_Excel/plot3.pdf")

plt.show()
plt.pause(0.1)

plt.figure()


hist_min = 0
hist_max = 60
hist_bins = np.linspace(hist_min,hist_max,3601)

hist_bins2 = []
for i in range(3601):
    hist_bins2.append(np.sqrt(i))



# counts,bins,bars = plt.hist(distOlin,bins=hist_bins)
# disthist = []
# for i in range(len(counts)):
#     counts[i] = counts[i]/((1+2*i)*np.pi)
#     for j in range(counts[i]):
#         disthist.append(i)
    
# counts2,bins,bars = plt.hist(distOlin2,bins=hist_bins)
# for i in range(len(counts)):
#     counts2[i] = counts2[i]/((1+2*i)*np.pi)



sns.set_style("white")
sns.distplot(distOlin, color="dodgerblue", label="Control", hist_kws={'alpha':.5}, kde_kws={'linewidth':2}, bins=hist_bins2)
sns.distplot(distOlin2, color="deeppink", label="Mutant", hist_kws={'alpha':.5}, kde_kws={'linewidth':2}, bins=hist_bins2)
# sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.yscale('log')

# plt.ylim(0,0.11)
plt.legend()

plt.savefig("E:/Workstation/Skeletonization_Excel/plot3.pdf")

plt.show()
plt.pause(0.1)











# plt.figure()
# plt.plot(bins[:-1],counts)
# plt.hist(bins[:-1],weights=counts,bins=bins[:-1])

# plt.figure()
# plt.plot(bins[:-1],counts2)
# plt.hist(bins[:-1],weights=counts2,bins=bins[:-1])

fig = plt.figure(figsize=(9.6,9.6))
ax = fig.add_subplot()


sns.set_style("white")
sns.distplot(bins[:-1],  color="dodgerblue", label="Control", hist_kws={'weights': counts,'alpha':.5}, kde_kws={'linewidth':0}, bins=hist_bins)
sns.distplot(bins[:-1],  color="deeppink", label="Mutant", hist_kws={'weights': counts2,'alpha':.0}, kde_kws={'linewidth':0}, bins=hist_bins)
# sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.ylim(0,0.11)
plt.legend()

# plt.savefig("E:/Workstation/Skeletonization_Excel/plot3.pdf")

plt.show()
plt.pause(0.1)


fig = plt.figure(figsize=(9.6,9.6))
ax = fig.add_subplot()


sns.set_style("white")
sns.distplot(bins[:-1],  color="dodgerblue", label="Control", hist_kws={'weights': counts,'alpha':.0}, kde_kws={'linewidth':0}, bins=hist_bins)
sns.distplot(bins[:-1],  color="deeppink", label="Mutant", hist_kws={'weights': counts2,'alpha':.5}, kde_kws={'linewidth':0}, bins=hist_bins)
# sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.ylim(0,0.11)
plt.legend()

# plt.savefig("E:/Workstation/Skeletonization_Excel/plot4.pdf")

plt.show()
plt.pause(0.1)



fig = plt.figure(figsize=(9.6,9.6))
ax = fig.add_subplot()

sns.set_style("white")
sns.distplot(bins[:-1],  color="dodgerblue", label="Control", hist_kws={'weights': counts,'alpha':.5}, kde_kws={'linewidth':0}, bins=hist_bins)
sns.distplot(bins[:-1],  color="deeppink", label="Mutant", hist_kws={'weights': counts2,'alpha':.5}, kde_kws={'linewidth':0}, bins=hist_bins)
# kwargs = dict(hist_kws={'alpha':.5}, kde_kws={'linewidth':0}, bins=hist_bins)
# sns.distplot(distOlin, color="dodgerblue", label="Control", **kwargs)
# sns.distplot(distOlin2, color="deeppink", label="Mutant", **kwargs)
# # sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.ylim(0,0.11)
plt.legend()

# plt.savefig("E:/Workstation/Skeletonization_Excel/plot5.pdf")

plt.show()
plt.pause(0.1)


fig = plt.figure(figsize=(9.6,9.6))
ax = fig.add_subplot()

hist_min = 0
hist_max = 60
hist_bins = np.linspace(hist_min,hist_max,61)

sns.set_style("white")
sns.distplot(bins[:-1],  color="dodgerblue", label="Control", hist_kws={'weights': counts,'alpha':.5}, kde_kws={'linewidth':2}, bins=hist_bins)
sns.distplot(bins[:-1],  color="deeppink", label="Mutant", hist_kws={'weights': counts2,'alpha':.5}, kde_kws={'linewidth':2}, bins=hist_bins)
# sns.distplot(x3, color="orange", label="minivan", **kwargs)
plt.xlim(0,60)
plt.ylim(0,0.11)
plt.legend()

# plt.savefig("E:/Workstation/Skeletonization_Excel/plot6.pdf")

plt.show()
plt.pause(0.1)


print("done")


def rotate(x,y,theta): #rotate x,y around xo,yo by theta (rad)
    xr=np.cos(theta)*(x)-np.sin(theta)*(y)
    yr=np.sin(theta)*(x)+np.cos(theta)*(y)
    return [xr,yr]
