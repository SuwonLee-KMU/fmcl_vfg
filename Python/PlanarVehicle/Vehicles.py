# Generated on: 2023-04-28 
# Author: Suwon Lee from Kookmin Univ.
# Vehicles class

import numpy as np
import PlanarVehicle as pv
from math import sin, cos

class Vehicle():
    def __init__(self, InitialState:pv.States = pv.States()):
        self.states = InitialState

    def dynamics(self, Input:pv.Inputs=pv.Inputs()):
        x, y, h, V = self.states.get_states()
        dx = V*cos(h)
        dy = V*sin(h)
        thrust, latax = Input.get_inputs()
        if V == 0:
            dh = 0
        else:
            dh = latax / V
        dV = thrust
        return np.array([dx, dy, dh, dV])
    
    
