# V_6 - Semi-Lagrangian Scheme
# To find the element of several coordinates
# and outside coordinates
# using neighbors elements

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

name_mesh = 'malha_cavity.msh'
number_equations = 3
mesh = trimsh.Linear('/home/marquesleandro/mesh',name_mesh, number_equations)
mesh.ien()
mesh.coord()

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""


#-----------------------------------------

x_d = np.zeros([mesh.npoints,1], dtype = float)
y_d = np.zeros([mesh.npoints,1], dtype = float)

vx = np.zeros([mesh.npoints,1], dtype = float)
vy = np.zeros([mesh.npoints,1], dtype = float)

dt = 0.005

for i in range(0,mesh.npoints):
 vx[i] = np.sin(mesh.y[i])
 vy[i] = 0.0

x_d = mesh.x
y_d = mesh.y

print '----------------'
print 'SEMI-LAGRANGIAN:'
print '----------------'

start_time = time()

nt = 300
for t in tqdm(range(0,nt)):
 save = InOut.Linear(x_d,y_d,mesh.IEN,mesh.npoints,mesh.nelem,vx,vx,vy,vx,vy)
 save.saveVTK('/home/marquesleandro/result','semi_lagrangian%s' %t) 

 for node in range(0,mesh.npoints):
  # mesh moving 
  x_d[node] = x_d[node] + vx[node]*dt
  y_d[node] = y_d[node] + vy[node]*dt

  x = float(x_d[node])
  y = float(y_d[node])

  length = []
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
  
    A = np.array([[x1,x2,x3],
                  [y1,y2,y3],
                  [1.0,1.0,1.0]])

    b = np.array([x,y,1.0])

    alpha = np.linalg.solve(A,b)
 
    if np.all(alpha >= 0.0) and np.all(alpha <= 1.0):
     ee = e + 175
     #print "elemento dominio %s" %ee
     #print "fazer interpolacao triangular"
     ww = 0
     break

    else:
     x_a = x1 - x
     x_b = x2 - x
     x_c = x3 - x
   
     y_a = y1 - y
     y_b = y2 - y
     y_c = y3 - y
 
     length1 = np.sqrt(x_a**2 + y_a**2)
     length2 = np.sqrt(x_b**2 + y_b**2)
     length3 = np.sqrt(x_c**2 + y_c**2)

     a_1 = [v1,length1]
     a_2 = [v2,length2]
     a_3 = [v3,length3]
 
     length.append(a_1)
     length.append(a_2)
     length.append(a_3)
   
     ww = 1
 
   # first neighbor is element found
   if ww == 0:
    break
 
   # coordinate doesn't found
   else:
    length_min = min(length, key=lambda k:k[1])
    node1 = node
    node = length_min[0]

    # outside domain
    if node == node1 and ww == 1:
     #node = node + 1
     #print "elemento contorno proximo ao no %s" %node
     #print "fazer regra da alavanca"    
     ww = 0
     break

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""

