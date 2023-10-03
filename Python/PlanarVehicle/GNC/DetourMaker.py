# Author: Suwon Lee from Kookmin University
# Date: April 29, 2023 1:54:24 PM GMT+9
# Code: POA 모듈을 사용해서 우회로(Detour)를 생성하는 클래스.

import POA
import POA.OffsetCurves as off
import numpy as np

class DetourMaker():
    def __init__(self, RefCurveObj:POA.Curves.PlanarPHCurve, num_sample_points=30):
        self.__Reference = RefCurveObj
        self.__sensed_obs = []
        self.num_sample_points = num_sample_points
        self.phi = off.ExponentialRBF(epsilon=self.num_sample_points)                                   # RBF 함수 설정 
        self.spt = off.SamplePoints(points=np.linspace(0,1,self.num_sample_points))                     # 샘플포인트 객체 생성     
        self.dlbc = off.DLBComputer_multipleObstacle(
            Reference=self.__Reference, SamplePointsObj=self.spt, step_length=0.05)
        self.spt.DLB = self.dlbc.determine_DLBs(self.__sensed_obs)                                      # 장애물과 샘플포인트 DLB 계산   
        self.oc = off.OffsetCurve(Reference=self.__Reference, SamplePoints=self.spt, RBF=self.phi)      # 오프셋곡선 객체 생성
        self.update_OffsetCurve()

    def __call__(self):
        pass

    def update_OffsetCurve(self):
        self.spt.DLB = self.dlbc.determine_DLBs(self.__sensed_obs)                          # 장애물과 샘플포인트 DLB 계산
        self.oc.SamplePoints = self.spt

    @property 
    def sensed_obs(self):
        return self.__sensed_obs
    @sensed_obs.setter 
    def sensed_obs(self, sensed_obs:list):
        self.__sensed_obs = sensed_obs
        self.update_OffsetCurve()