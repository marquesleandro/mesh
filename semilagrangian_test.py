# V_2 - Semi-Lagrangian Scheme
# To find the element of several coordinates
# and outside coordinates with convex surface

# =======================
# Importing the libraries
# =======================

import sys
sys.path.insert(0, '../lib_class')

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg
import scipy.linalg
import import_msh
from tqdm import tqdm
from time import time


print '------------'
print 'IMPORT MESH:'
print '------------'

start_time = time()

name_mesh = 'semilagrangian_test2.msh'
number_equations = 1
mesh = import_msh.Linear2D('/home/marquesleandro/mesh',name_mesh, number_equations)
mesh.coord()
mesh.ien()

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""

#-----------------------------------------

print '----------------'
print 'SEMI-LAGRANGIAN:'
print '----------------'

start_time = time()



# varios pontos
points1_x = np.zeros([6,1], dtype = float)
points1_y = np.zeros([6,1], dtype = float)

#no 0 -> proximo ao 11
points1_x[0] = 1.1
points1_y[0] = 0.4


#no 1 -> proximo ao 23
points1_x[1] = 0.4
points1_y[1] = 0.3

#no 2 -> proximo ao 5
points1_x[2] = 0.0
points1_y[2] = 1.0

#no 3 -> proximo ao 32
points1_x[3] = 0.8
points1_y[3] = 0.1

#no 4 -> proximo ao 12
points1_x[4] = 1.1
points1_y[4] = 0.2

#no 5 -> proximo ao 28
points1_x[5] = 0.1
points1_y[5] = 0.9

for i in range(0,6): #range mesh.npoints
 x = float(points1_x[i])
 y = float(points1_y[i])
 
 node = i
 length = []
 barycentric = []
 element = []
 breaking = 0

 while breaking == 0:
  for j in mesh.far_neighbors_nodes[node]:
   xx = float(mesh.x[j])
   yy = float(mesh.y[j])
 
   x_a = xx - x
   y_a = yy - y
  
   length1 = np.sqrt(x_a**2 + y_a**2)

   a_1 = [j,length1]
 
   length.append(a_1)
   
  length_min = min(length, key=lambda k:k[1])
  node1 = node
  node = length_min[0]

  # node more short
  if node == node1:
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

    barycentric.append(alpha)
    element.append(e)

   barycentric = np.array(barycentric)
   idx = (barycentric >= 0.0) & (barycentric <= 1.0)
   aa = np.where((idx == True).all(axis=1))
   print idx
   print aa

   # inside element
   try:
    aa = int(aa[0])
    e = element[aa]

    v1 = mesh.IEN[e][0]
    v2 = mesh.IEN[e][1]
    v3 = mesh.IEN[e][2]

    x1 = float(mesh.x[v1])
    x2 = float(mesh.x[v2])
    x3 = float(mesh.x[v3])

    y1 = float(mesh.y[v1])
    y2 = float(mesh.y[v2])
    y3 = float(mesh.y[v3])
  

    A1 = 0.5*np.linalg.det(np.array([[1, x, y],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
 
    A2 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x, y],
                                     [1, x3, y3]]))
 
    A3 = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x, y]]))
 
    At = 0.5*np.linalg.det(np.array([[1, x1, y1],
                                     [1, x2, y2],
                                     [1, x3, y3]]))
   
    L1 = A1/At
    L2 = A2/At
    L3 = A3/At
     
    N1 = L1
    N2 = L2
    N3 = L3

    
    print "interpolacao no %s" %node
    print ""

    breaking = 1
    break

    # outside domain
   except TypeError:
    print "fora dominio no %s" %node
    print ""
    
    breaking = 1
    break

  # node far yet
  else:
   breaking = 0



end_time = time()
print ""
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""


