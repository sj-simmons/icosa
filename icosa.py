import pygame, sys
from pygame.locals import *
import numpy as np
from math import sin, cos, pi
from functools import reduce

winWidth = 800    
winHeight = 600 
FPS = 120 

WHITE = (255, 255, 255)
GREY = (180, 180, 180)
BLACK = (0, 0, 0)

''' A class for working with homogeneous coordinates of 3d points'''
class Pmatrix(np.ndarray):
  def __new__(cls, input_array):
    obj = np.asarray(input_array).view(cls)
    return obj

  def homog(self):
    """ return the matrix with the points defined by self's columns, 
        homogenized""" 
    return Pmatrix(np.vstack((self,np.ones((1,self.shape[1])))))

  def dehom(self):
    """ return the matrix with points defined by selfs's columns, 
        dehomogenized""" 
    A = self[:self.shape[0]-1,:]
    for j in range(A.shape[1]):
      for i in range(A.shape[0]):
        A[i,j] /= self[self.shape[0]-1,j]
    return A 

  ''' Project onto the xy-plane along the z axis from the point (0,0,d)'''
  def Pxy(d):
    return Pmatrix(np.array( [[ 1, 0,    0, 0],
                              [ 0, 1,    0, 0],
                              [ 0, 0,    0, 0],
                              [ 0, 0, -1/d, 1]]))

  def Rxy(theta):
    return Pmatrix( np.array( [[ cos(theta*pi/180), -sin(theta*pi/180), 0, 0],\
                               [ sin(theta*pi/180), cos(theta*pi/180), 0, 0],\
                               [ 0, 0, 1, 0],\
                               [ 0, 0, 0, 1]]))

  def Ryz(theta):
    return Pmatrix( np.array( [[ 1, 0, 0, 0],\
                               [ 0, cos(theta*pi/180), -sin(theta*pi/180), 0],\
                               [ 0, sin(theta*pi/180),  cos(theta*pi/180), 0],\
                               [ 0, 0, 0, 1]]))

  def Rxz(theta):
    return Pmatrix( np.array( [[ cos(theta*pi/180), 0, -sin(theta*pi/180), 0],\
                               [ 0, 1, 0, 0],\
                               [ sin(theta*pi/180), 0,  cos(theta*pi/180), 0],\
                               [ 0, 0, 0, 1 ]]))

  def T(x,y,z):
    return Pmatrix(np.array( [[ 1, 0, 0, x],
                              [ 0, 1, 0, y],
                              [ 0, 0, 1, z],
                              [ 0, 0, 0, 1]]))

  """ return the matrix product of the arguments """
  def xdot(*args):
    return reduce(Pmatrix.dot, args)

  def twoDverts(self):
    """ return columns of 3 x m self as tuples in R^2 by stripping off last 
        coordinate"""
    V = self[:self.shape[0]-1,:]
    tup_lst = []
    for j in range(V.shape[1]):
      tup_lst.append(tuple(V[:,j]))
    return tup_lst 

  def verts2Pmatrix(verts):
    """ return a Pmatr in which each column is a point version of the tuple in 
        verts, in homogeneous coords"""
    return Pmatrix(np.array(verts).reshape(len(verts),3).transpose())

  def min_3d_dists(verts):
    """ return 2-tuples corresponding to edges with minimal length"""
    edges = []
    dists = {} 
    for i in range(len(verts)):
      for j in range(len(verts)):
        if j > i:
          dists[i,j] = (verts[i][0]-verts[j][0])**2+(verts[i][1]-verts[j][1])**\
                                               2+(verts[i][2]-verts[j][2])**2
    min_dist = min(dists.values())
    for key in dists.keys():
      if dists[key] - min_dist < .1:
        edges.append(key)
    return edges

class Polytope():
  def __init__(self,verts, edges):
    self.verts = verts
    self.edges = edges

class Icosa(Polytope):
  r = 1.61803398875  # golden ratio
  l = 70
  # verts of a regular icosahedron
  rect1 = [(r*l,l,0), (r*l,-l,0), (-r*l,-l,0), (-r*l,l,0)] #xy-plane\   scaled
  rect2 = [(l,0,r*l), (l,0,-r*l), (-l,0,-r*l), (-l,0,r*l)] #xz       |- golden 
  rect3 = [(0,r*l,l), (0,r*l,-l), (0,-r*l,-l), (0,-r*l,l)] #yz      /   rects
  verts = rect1 + rect2 + rect3

  edges = Pmatrix.min_3d_dists(verts)

  V = Pmatrix.verts2Pmatrix(verts).homog()

  def __init__(self, size = 70):
    Icosa.l = size
    super().__init__(Icosa.verts, Icosa.edges)

  def draw(self):
    d = 700  
    P = Pmatrix.Pxy(d)
    R = np.dot(Pmatrix.Ryz(2/30),Pmatrix.Rxz(8/30))
    Txy = Pmatrix.T(winWidth/3,winHeight/3,0)
    Tz = Pmatrix.T(0,0,300)
    Icosa.V = np.dot(R,Icosa.V)
    if shadow:
      VV = Pmatrix.xdot(P,Txy,Tz,Icosa.V).dehom()
      verts = VV.twoDverts()
      for edge in Icosa.edges:
        pygame.draw.line(screen,GREY,verts[edge[0]],verts[edge[1]],6)
    V = Pmatrix.xdot(Txy,Tz,P,Icosa.V).dehom()
    verts = V.twoDverts()
    for edge in Icosa.edges:
      pygame.draw.line(screen,BLACK,verts[edge[0]],verts[edge[1]],4)
      
pygame.init()
fpsClock = pygame.time.Clock()
screen = pygame.display.set_mode((winWidth,winHeight))
pygame.display.set_caption('Icosahedron demo')

icosa = Icosa(300)
shadow = False

while True:
  for event in pygame.event.get():
    if event.type == QUIT or (event.type == KEYUP and event.key == K_q):
      pygame.quit()
      sys.exit()
    if event.type == KEYUP and event.key == K_s:
      shadow = not shadow
    if event.type == KEYUP and event.key == K_p:
      print(Icosa.V) 

  screen.fill(WHITE)
  
  icosa.draw()

  pygame.display.update()
  fpsClock.tick(FPS)
