import sys
from pyMCDS_optional_meshes import pyMCDS
import numpy as np

print(sys.argv)
frame = int(sys.argv[1])

fname = "output%08d.xml" % frame
mcds = pyMCDS(fname,'.')  
tmins = mcds.get_time()
print('time (mins)=',tmins)
print('time (days)=',tmins/1440.)

print("[discrete_cells].keys() = ",mcds.data['discrete_cells'].keys())
# print("IDs= ",mcds.data['discrete_cells']['ID'])
ncells = len(mcds.data['discrete_cells']['ID'])
print('num cells = ',ncells)

#xyz = np.empty((ncells,3))
xyz = np.zeros((ncells,3))
xvals = mcds.data['discrete_cells']['position_x']
#print("position_x= ",mcds.data['discrete_cells']['position_x'])
yvals = mcds.data['discrete_cells']['position_y']
zvals = mcds.data['discrete_cells']['position_z']
print("x position range: ",xvals.min(),xvals.max())
print("y position range: ",yvals.min(),yvals.max())
print("z position range: ",zvals.min(),zvals.max())

xorien = mcds.data['discrete_cells']['orientation_x']
#print("position_x= ",mcds.data['discrete_cells']['position_x'])
yorien = mcds.data['discrete_cells']['orientation_y']
zorien = mcds.data['discrete_cells']['orientation_z']
print("x orientation range: ",xorien.min(),xorien.max())
print("y orientation range: ",yorien.min(),yorien.max())
print("z orientation range: ",zorien.min(),zorien.max())

axis_a = mcds.data['discrete_cells']['axis_a']
axis_b = mcds.data['discrete_cells']['axis_b']
axis_c = mcds.data['discrete_cells']['axis_c']
print("axis_a range: ",axis_a.min(),axis_a.max())
print("axis_b range: ",axis_b.min(),axis_b.max())
print("axis_c range: ",axis_c.min(),axis_c.max())

# "_n" where n=[0,T) (# of cell types)
# attack_rates= mcds.data['discrete_cells']['attack_rates']
# print("attack_rates_0 =",attack_rates)
