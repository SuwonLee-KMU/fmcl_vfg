# Generated on: 2023-04-28 
# Author: Suwon Lee from Kookmin Univ.

import numpy as np
from math import cos, sin

class States():
    def __init__(self):
        self.x = 0
        self.y = 0
        self.heading = 0
        self.speed = 0

    def get_states(self):
        return np.array([self.x, self.y, self.heading, self.speed])
    
    def set_states(self, X):
        self.x, self.y, self.heading, self.speed = X

    def get_velocity(self):
        u = self.speed*cos(self.heading)
        v = self.speed*sin(self.heading)
        return u, v


    
