# =======================
# Importing the libraries
# =======================

import sys
sys.path.insert(0, '../lib_class')

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg
import scipy.linalg
import trimsh
import trielem
from tricond import b_bc
import InOut
from tqdm import tqdm
from time import time


print '------------'
print 'IMPORT MESH:'
print '------------'

start_time = time()

name_mesh = 'semilagrangian_test.msh'
number_equations = 3
mesh = trimsh.Linear('/home/marquesleandro/fem/mesh',name_mesh, number_equations)
mesh.ien()
mesh.coord()

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""


points1_x = np.zeros([3,1], dtype = float)
points1_y = np.zeros([3,1], dtype = float)

#200
points1_x[0] = 0.5
points1_y[0] = 0.4

#187
points1_x[1] = 2.2
points1_y[1] = 0.7

#261
points1_x[2] = 3.0
points1_y[2] = 0.5

# node 62
node = 62
x = 1.9
y = 0.5


element_avg = []

barycenter = np.array([10.0,10.0,10.0])
ww = 1

while ww == 1:
 for e in mesh.neighbors_elements[node]:
  v1 = mesh.IEN[e][0]
  v2 = mesh.IEN[e][1]
  v3 = mesh.IEN[e][2]

  x1 = float(mesh.x[v1])
  x2 = float(mesh.x[v2])
  x3 = float(mesh.x[v3])

  y1 = float(mesh.y[v1])
  y2 = float(mesh.y[v2])
  y3 = float(mesh.y[v3])
  
#  x = float(points1_x[i])
#  y = float(points1_y[i])

  A = np.array([[x1,x2,x3],
                [y1,y2,y3],
                [1.0,1.0,1.0]])

  b = np.array([x,y,1.0])

  alpha = np.linalg.solve(A,b)

  barycenter = ((barycenter,alpha))

  avg = np.average(np.sqrt(alpha**2))
  aa2 = [e, avg]
  element_avg.append(aa2) 

 near_element = min(element_avg, key=lambda k:k[1])
 near_element = near_element[0]
  
 v1 = mesh.IEN[near_element][0]
 v2 = mesh.IEN[near_element][1]
 v3 = mesh.IEN[near_element][2]

 x1 = float(mesh.x[v1])
 x2 = float(mesh.x[v2])
 x3 = float(mesh.x[v3])

 y1 = float(mesh.y[v1])
 y2 = float(mesh.y[v2])
 y3 = float(mesh.y[v3])
  
 x_a = x1 - x
 x_b = x2 - x
 x_c = x3 - x
  
 y_a = y1 - x
 y_b = y2 - x
 y_c = y3 - x
 
 length1 = np.sqrt(x_a**2 + y_a**2)
 length2 = np.sqrt(x_b**2 + y_b**2)
 length3 = np.sqrt(x_c**2 + y_c**2)

 lenght = [[v1,length1],[v2,length2],[v3,length3]]
 lenght_min = min(lenght, key=lambda k:k[1])
 node = lenght_min[0]
 print node
 print element_avg

 for i in range(0,len(barycenter)):
  if np.all(barycenter[i] >= 0.0) and np.all(barycenter[i] <= 1.0):
   ee = e + 175
   print ee
   ww = 0
 
  else:
   ww = 1

 element_avg = []
 barycenter = np.array([10.0,10.0,10.0])
 print element_avg
