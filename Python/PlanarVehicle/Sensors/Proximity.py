# Author: Suwon Lee from Kookmin University
# Date: April 29, 2023 11:36:12 AM GMT+9
# Code: Proximity Sensor Class for PlanarVehicle Module
# 모든 장애물의 위치를 사전에 받아서 가지고 있으며,
# 장애물이 반경내로 들어오는 경우 이를 알려준다.
# 실제 물리적인 센서를 모델링하는 것은 아님.

import PlanarVehicle as pv
import numpy as np
from math import sqrt

class Proximity():
    def __init__(self, 
                 radius:float=0.01, 
                 States:pv.States=pv.States(),
                 obs_list:list=[],
                 ):
        """
        본 클래스의 인스턴스 속성값에 대한 설명은 아래와 같다.
        radius : 센서의 탐지반경
        States : 센서가 탑재된 PlanarVehicle의 상태변수 객체
        obs_list : 전체 임무영역의 모든 장애물들의 리스트.
        sensing_list : 각 장애물마다에 대한 감지여부 리스트.
        """
        self.__radius = radius
        self.__States = States
        self.__obs_list = obs_list # list of POA.Obstacles.CircularObstacles 
        self.__sensing_list = [False]*len(self.obs_list)
        self.is_updated = False

    def update_sensing_list(self):
        """
        현재 상태변수에 대한 sensing_list 속성값을 업데이트한다.
        """
        px, py = self.__States.x, self.__States.y
        is_sensed = []
        for id, obs in enumerate(self.obs_list):
            cx, cy = obs.center
            center_dist = sqrt((px-cx)**2 + (py-cy)**2)
            if center_dist < (self.__radius + obs.radius):
                is_sensed.append(True)
            else:
                is_sensed.append(False)
        if self.__sensing_list == is_sensed:
            self.is_updated = False
        else:
            self.is_updated = True
            self.__sensing_list = is_sensed

    @property
    def radius(self):
        return self.__radius
    @radius.setter
    def radius(self, radius:float):
        self.__radius = radius
    @property
    def States(self):
        return self.__States
    @States.setter
    def States(self, States:float):
        self.__States = States
    @property
    def obs_list(self):
        return self.__obs_list
    @obs_list.setter
    def obs_list(self, obs_list:float):
        self.__obs_list = obs_list
    @property
    def sensing_list(self):
        return self.__sensing_list
    @sensing_list.setter
    def sensing_list(self, sensing_list:float):
        self.__sensing_list = sensing_list

    