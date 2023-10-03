# Generated on: 2023-04-28 
# Author: Suwon Lee from Kookmin Univ.

import PlanarVehicle as pv
from math import cos, sin

class Simout():
    def __init__(self):
        self.x = []
        self.y = []
        self.heading = []
        self.speed = []
        self.thrust = []
        self.latax = []
        self.u = []
        self.v = []

    def add_states(self, states:pv.States, inputs:pv.Inputs):
        x, y, h, V = states.get_states()
        t, l = inputs.get_inputs()
        self.x.append(x)
        self.y.append(y)
        self.heading.append(h)
        self.speed.append(V)
        self.thrust.append(t)
        self.latax.append(l)
        self.u.append(V*cos(h))
        self.v.append(V*sin(h))

        

class Simout_VFG():
    def __init__(self):
        self.u = []
        self.v = []

    def add_states(self, u, v):
        self.u.append(u)
        self.v.append(v)


class Simout_Proximity():
    def __init__(self, obs_list:list=[]):
        self.obs_list = obs_list
        self.stacked_sensing_list = []

    def add_states(self, sensing_list):
        self.stacked_sensing_list.append(sensing_list)