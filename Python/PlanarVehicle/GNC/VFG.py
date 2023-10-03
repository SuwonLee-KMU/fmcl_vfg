# Generated on: 2023-04-28 
# Author: Suwon Lee from Kookmin Univ.

import PlanarVehicle as pv
import numpy as np
from math import atan

class VectorFieldGuidance():
    def __init__(self, CurveObject, convergence_factor, gain, States:pv.States=pv.States()):
        """
        인스턴스 속성값에 대한 설명
        CurveObject : .__call__(), .derivative() 메서드를 구현해야 한다.
        """
        self.CurveObject = CurveObject
        self.convergence_factor = convergence_factor
        self.gain = gain
        self.__States = States
        self.vel_cmd = None
        self.latax_cmd = None
        
    def __call__(self):
        pass

    def compute_vel_cmd(self):
        x, y, h, V = self.__States.get_states()
        pos_tgt, tau_tgt, distance = self.nearest_point_from(x, y)
        tan, nor = self.basis_vectors_at(tau_tgt)
        vec_conv = (pos_tgt - np.array([x,y])) / np.linalg.norm(pos_tgt - np.array([x,y]))
        vec_trav = tan
        vel_vector = np.reshape(0.3*vec_conv + 0.6*vec_trav, newshape=[2])
        return vel_vector/np.linalg.norm(vel_vector)

    def compute_latax_cmd(self):
        vel = self.__States.get_velocity()
        error = np.arctan2(self.vel_cmd[1]*vel[0]-self.vel_cmd[0]*vel[1],self.vel_cmd[0]*vel[0]+self.vel_cmd[1]*vel[1]) # angle between 2D vectors
        latax_cmd = error*self.gain
        return latax_cmd

    @property 
    def States(self):
        return self.__States
    @States.setter 
    def States(self, States:pv.States):
        self.__States = States
        self.vel_cmd = self.compute_vel_cmd()
        self.latax_cmd = self.compute_latax_cmd()

    def basis_vectors_at(self, tau):
        drv = np.transpose(np.reshape(self.CurveObject.derivative(tau),newshape=[2,-1]))
        tan = np.array([d/np.linalg.norm(d+0.0000001) for d in drv])
        nor = np.array([[-t[1],t[0]] for t in tan])
        return tan, nor
    
    def dist_based_coef(self, dist):
        k_conv = atan(dist*self.convergence_factor)/np.pi*2
        k_trav = np.sqrt(1-k_conv**2)
        coef   = [k_conv, k_trav]
        return coef

    def nearest_point_from(self, x, y):
        tau = np.linspace(0,1,201)
        pos = np.transpose(self.CurveObject(tau))
        norms = [np.linalg.norm(np.reshape(p,newshape=[2])-np.array([x,y])) for p in pos]
        minval = min(norms)
        index = norms.index(minval)
        return pos[index], tau[index], minval
    
    
