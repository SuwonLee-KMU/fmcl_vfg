% Author: Suwon Lee from Kookmin University
% Date: February 7, 2023 5:18:52 PM GMT+9
% Code: VFG using Polynomial curve as a reference path.

clear all; clc;

%% 1. 3차원 경로 데이터 불러오기.
data = load("data0.mat");
cvs  = data.curves;
for i = 1:numel(cvs)
    scp{i} = SpaceCurvePoly(cvs{i}.coef_x, cvs{i}.coef_y, cvs{i}.coef_z);
end
mcp = MultipleCurvesPoly(scp); % 다수 곡선을 하나의 매개변수 tau(0~1)로 연결해주는 클래스 

%% 2. 불러온 곡선의 위치를 계산 (시각화 목적)
tau = linspace(0,1,201);    % 곡선의 매개변수 정의.
pos = mcp.feval(tau);       % 곡선의 위치. [ndata x 3]
tan = mcp.tangent(tau);     % 곡선의 접선벡터. [ndata x 3]

%% 3. 벡터필드 유도기법 
vfg  = VFGParametric(mcp, 10);                              % 벡터필드 유도기법 객체 생성
rpos = vfg.generate_random_position_near_curve(200,0.1);    % 곡선 근처의 무작위 점 생성
v_d  = vfg.feval(rpos);                                     % 생성된 점에서 곡선에 수렴하는 속도벡터 명령 계산

%% 4. 시각화
fig = figure(1); clf;
plot3(pos(:,1),pos(:,2),pos(:,3),...
    'linewidth',3, 'displayname',"Ref. path");
hold on;
quiver3(rpos(:,1),rpos(:,2),rpos(:,3),v_d(:,1),v_d(:,2),v_d(:,3),...
    "displayname", "Guidance command(v_d)");
scatter3(rpos(:,1),rpos(:,2),rpos(:,3),...
    'markerfacecolor',[1,0.7,0.7],...
    'markeredgecolor','none',...
    "displayname","Random points");
legend("location", "northeast")
grid on; box on; axis equal;