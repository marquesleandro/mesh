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

name_mesh = 'semilagrangian_test2.msh'
number_equations = 1
mesh = trimsh.Linear('/home/marquesleandro/mesh',name_mesh, number_equations)
mesh.ien()
mesh.coord()

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""

#-----------------------------------------

print '----------------'
print 'SEMI-LAGRANGIAN:'
print '----------------'

start_time = time()



# varios pontos
points1_x = np.zeros([9,1], dtype = float)
points1_y = np.zeros([9,1], dtype = float)

#node 0
#238 elemento vizinho proximo do ponto 1
points1_x[0] = 0.23
points1_y[0] = 0.97

#node 1
#490 elemento vizinho distante do ponto 2
points1_x[1] = 0.24
points1_y[1] = 0.44

#node 2
#334 elemento mega distante do ponto 3
points1_x[2] = 0.14
points1_y[2] = 0.12

#node 3
#277 elemento vizinho mega distante do ponto 4 colado no contorno
points1_x[3] = 0.0
points1_y[3] = 0.21

#node 4
#noh 5 no vertice imovel - elemento 587
points1_x[4] = 0.0
points1_y[4] = 0.0

#node 5
#noh vizinho proximo do ponto 6 fora do dominio
points1_x[5] = 0.02
points1_y[5] = 1.05

#node 6
#noh 7 no contorno imovel - elemento 51
points1_x[6] = 0.25
points1_y[6] = 0.95

#node 7
#285 elemento mega distante do ponto 8 com uma superficie convexa
points1_x[7] = 0.67
points1_y[7] = 0.47

#node 8
#noh mega distante proximo do ponto 36 fora do dominio com uma superficie convexa
points1_x[8] = 1.2
points1_y[8] = 0.15



for node in range(0,9): #range mesh.npoints
 x = float(points1_x[node])
 y = float(points1_y[node])

 length = []
 ww = 1
 print ""
 print node

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
    ee = e + 11
    print "elemento dominio %s" %ee
    print "fazer interpolacao triangular" 
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
   print node

   # outside domain
   if node == node1 and ww == 1:
    node = node + 1
    print "elemento contorno proximo ao no %s" %node
    print "fazer regra da alavanca"
    ww = 0
    break


end_time = time()
print ""
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""


