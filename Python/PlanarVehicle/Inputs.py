# Generated on: 2023-04-28 
# Author: Suwon Lee from Kookmin Univ.

import numpy as np

class Inputs():
    def __init__(self, ):
        self.thrust = 0
        self.latax = 0

    def get_inputs(self):
        return np.array([self.thrust, self.latax])
    
    def set_inputs(self, u):
        self.thrust, self.latax = u


    