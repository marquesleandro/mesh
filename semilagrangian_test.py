# V_5 - Semi-Lagrangian Scheme
# To find the element of several coordinates
# and outside coordinates

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
mesh = trimsh.Linear('/home/marquesleandro/mesh',name_mesh, number_equations)
mesh.ien()
mesh.coord()

end_time = time()
print 'time duration: %.1f seconds' %(end_time - start_time)
print ""

#-----------------------------------------
# varios pontos
points1_x = np.zeros([9,1], dtype = float)
points1_y = np.zeros([9,1], dtype = float)

#node 0
#204 elemento vizinho proximo do ponto 1
points1_x[0] = 0.1
points1_y[0] = 0.8

#node 1
#180 elemento vizinho distante do ponto 2
points1_x[1] = 2.0
points1_y[1] = 0.7

#node 2
#3 noh no contorno imovel
points1_x[2] = 8.275
points1_y[2] = 1.0

#node 3
#4 noh no vertice imovel
points1_x[3] = 10.0
points1_y[3] = 1.0

#node 4
#5 noh vizinho proximo do ponto 5 fora do dominio
points1_x[4] = 10.1
points1_y[4] = 0.1

#node 5
#7 noh vizinho distante do ponto 6 fora do dominio
points1_x[5] = 0.5
points1_y[5] = 1.1

#node 6
#36 noh vizinho mega distante do ponto 7 fora do dominio
points1_x[6] = 8.4
points1_y[6] = -0.1

#node 7
#228 elemento vizinho mega distante do ponto 8
points1_x[7] = 9.2
points1_y[7] = 0.7

#node 8
#317 elemento vizinho mega distante do ponto 9 colado no contorno
points1_x[8] = 5.8
points1_y[8] = 0.1



'''
#node
node = 61
x = 9.1
y = 0.5
'''


for node in range(0,9): #range mesh.npoints
 x = float(points1_x[node])
 y = float(points1_y[node])

 length = []
 ww = 1
 
 print node

#length = []
#ww = 1
 while ww == 1:
#for i in range(0,3):

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
 
  if ww == 0:
    break
  
  else:
   length_min = min(length, key=lambda k:k[1])
   node1 = node
   node = length_min[0]

   if node == node1 and ww == 1:
    node = node + 1
    print "elemento contorno proximo ao no %s" %node
    print "fazer regra da alavanca"
    ww = 0
    break

