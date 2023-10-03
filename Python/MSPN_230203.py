# Generated on: 2023-01-31 
# Author: Suwon Lee from Kookmin Univ.

import numpy as np
import sympy
import math
# import plotly.graph_objects as go
# import scipy

def normalize(vector):
    if len(np.shape(vector)) == 1:
        return vector / np.linalg.norm(vector)
    else:
        raise ValueError(f'Invalid shape of the input vector : {np.shape(vector)}')

class HermiteInterpolation():
    # PH Quintic을 사용해서 평면상의 Hermite Interpolation을 수행하는 클래스
    def __init__(self, Ps, Ts, scaling, tau=None, mode=0):
        Pi, Pf = Ps
        Ti, Tf = Ts
        if tau is None:
            ndata = 100
            tau = np.linspace(0,1,ndata)
        Ti = np.reshape(scaling*normalize(np.reshape(Ti,newshape=[2,])),newshape=[2])
        Tf = np.reshape(scaling*normalize(np.reshape(Tf,newshape=[2,])),newshape=[2])
        Delta_p0 = Ti
        Delta_p4 = Tf
        p1 = Pi + Ti
        p4 = Pf - Tf
        Delta_p4p0 = p4 - p1
        if mode == 0:
            sign1 = 1
            sign2 = 1
        elif mode == 1:
            sign1 = -1
            sign2 = 1
        elif mode == 2:
            sign1 = 1
            sign2 = -1
        else:
            sign1 = -1
            sign2 = -1
        self.sign1 = sign1
        self.sign2 = sign2
            
        self.us, self.vs    = HermiteInterpolation.get_uvs(Delta_p0, Delta_p4, Delta_p4p0, sign1, sign2)
        self.cps            = HermiteInterpolation.get_control_points(self.us, self.vs, Pi)
        self.sigma_coef     = self.get_sigma_coefficients()
        self.parametric_speed = HermiteInterpolation.Bernstein_polynomial(self.sigma_coef)
        self.Pi = Pi
        self.Pf = Pf
        self.Ti = Ti
        self.Tf = Tf
        self.tau = tau

        # Bezier curve
        betas_x = [self.cps[i,0] for i in range(0,6)]
        betas_y = [self.cps[i,1] for i in range(0,6)]
        self.ph_x = HermiteInterpolation.Bernstein_polynomial(betas_x)
        self.ph_y = HermiteInterpolation.Bernstein_polynomial(betas_y)

    def visualize(self, fig=None):
        if fig is None:
            fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=self.ph_x(self.tau),
            y=self.ph_y(self.tau),
            name='PH curve'
        ))
        betas_x = [self.cps[i,0] for i in range(0,6)]
        betas_y = [self.cps[i,1] for i in range(0,6)]
        fig.add_trace(go.Scatter(
            x = betas_x,
            y = betas_y,
            line = dict(dash='dash'),
            name='Control points'
        ))
        return fig
        
    # Bernstein basis polynomial
    @staticmethod
    def Bernstein_basis_polynomial(n, nu):
        def b_n_nu(x):
            return math.comb(n, nu) * x**nu * (1-x)**(n-nu)
        return b_n_nu

    # Bernstein Polynomial
    @staticmethod
    def Bernstein_polynomial(betas):
        n = np.size(betas) - 1 # n: degree of the polynomial
        def Bn(x):
            basis = np.array([HermiteInterpolation.Bernstein_basis_polynomial(n,nu)(x) for nu in range(n+1)])
            return np.dot(betas, basis)
        return Bn

    @staticmethod
    def nonzero_sign(x):
        return 2*(x>=0).astype(int)-1

    # 식 (25.4)
    @staticmethod
    def get_uvs(Delta_p0, Delta_p4, Delta_p4p0, sign1, sign2):
        dx0, dy0 = Delta_p0
        dx4, dy4 = Delta_p4
        u0 = sign1*np.sqrt(5/2)*np.sqrt(np.linalg.norm(Delta_p0)+dx0)
        v0 = sign1*np.sqrt(5/2)*HermiteInterpolation.nonzero_sign(dy0)*np.sqrt(np.linalg.norm(Delta_p0)-dx0)
        u2 = sign2*np.sqrt(5/2)*np.sqrt(np.linalg.norm(Delta_p4)+dx4)
        v2 = sign2*np.sqrt(5/2)*HermiteInterpolation.nonzero_sign(dy4)*np.sqrt(np.linalg.norm(Delta_p4)-dx4)
        a = 9/16*(u0**2-v0**2+u2**2-v2**2)+5/8*(u0*u2-v0*v2)+15/2*(Delta_p4p0[0])
        b = 9/8*(u0*v0+u2*v2)+5/8*(u0*v2+u2*v0)+15/2*(Delta_p4p0[1])
        c = np.sqrt(a**2+b**2)
        u1 = -3/4*(u0+u2) + 1/np.sqrt(2)*np.sqrt(c+a)
        v1 = -3/4*(v0+v2) + np.sign(b)*1/np.sqrt(2)*np.sqrt(c-a)
        return [u0, u1, u2], [v0, v1, v2]

    # 식 (17.6)
    @staticmethod
    def get_control_points(us, vs, p0):
        p0 = np.reshape(p0, newshape=[2])
        u0, u1, u2 = us
        v0, v1, v2 = vs
        p1 = p0 + 1/5*np.array([u0**2-v0**2, 2*u0*v0])
        p2 = p1 + 1/5*np.array([u0*u1-v0*v1, u0*v1+u1*v0])
        p3 = p2 + 2/15*np.array([u1**2-v1**2, 2*u1*v1]) + 1/15*np.array([u0*u2-v0*v2,u0*v2+u2*v0])
        p4 = p3 + 1/5*np.array([u1*u2-v1*v2, u1*v2+u2*v1])
        p5 = p4 + 1/5*np.array([u2**2-v2**2, 2*u2*v2])
        return np.array([p0, p1, p2, p3, p4, p5])

    # 식 (17.12)
    def get_sigma_coefficients(self):
        u0, u1, u2 = self.us
        v0, v1, v2 = self.vs
        sigma       = [None]*5
        sigma[0]    = u0**2 + v0**2
        sigma[1]    = u0*u1 + v0*v1
        sigma[2]    = 2/3*(u1**2 + v1**2) + 1/3*(u0*u2 + v0*v2)
        sigma[3]    = u1*u2 + v1*v2
        sigma[4]    = u2**2 + v2**2
        return sigma

    def eval_sigma(self, tau=None):
        if tau is None:
            tau = self.tau
        return np.array(self.parametric_speed(tau))

    def eval(self, tau=None):
        if tau is None:
            tau = self.tau
        return np.array([self.ph_x(tau), self.ph_y(tau)])


class PreimageComplexPolynomial():
    def __init__(self, coef):
        self.coefficients = coef # Array of complex numbers. 낮은 차수(0차)부터 쓴다.

    # 입력한 계수에 대해 선형 복소방정식 함수를 출력하는 함수.
    def get_linear_complex_polynomial_handle(self):
        def func(z=0):
            dim = len(list(self.coefficients))
            zbasis = np.array([z**i for i in range(dim)])
            return np.dot(self.coefficients, zbasis)
        return func

    def __call__(self, z):
        func = self.get_linear_complex_polynomial_handle()
        return func(z)

    def symbolic(self):
        z = sympy.Symbol('z')
        func = self.get_linear_complex_polynomial_handle()
        return sympy.collect(func(z),z)


class MSPN(): # Minimum Surface with Pythagorean Normal
    def __init__(self, preimage1, preimage2, preimage3, r0 = [0,0,0]):
        self.u  = preimage1
        self.v  = preimage2
        self.w  = preimage3
        self.r0 = np.array(r0)

    def __call__(self, z):
        func = self.get_surface_handle(form='complex')
        return func(z)

    # u, v, w 선형복소함수에 대해 Enneper–Weierstrass parameterization을 수행하여, 복소함수를 출력하는 함수.
    def get_holomorphic_handles(self):
        u, v, w = self.u, self.v, self.w
        def Phi1(z):
            return w(z)*(u(z)**2+v(z)**2)
        def Phi2(z):
            return 2*w(z)*u(z)*v(z)
        def Phi3(z):
            return 1j*w(z)*(u(z)**2+v(z)**2)
        return Phi1, Phi2, Phi3

    # 다항식의 계수 u, v, w에 대해 Enneper–Weierstrass parameterization을 수행하여, 복소함수의 계수를 출력하는 함수.
    # 계수 리스트는 높은차수부터 낮은차수 순서로 쓴다.
    def get_holomorphic_coefs(self):
        u = np.flip(self.u.coefficients,axis=0)
        v = np.flip(self.v.coefficients,axis=0)
        w = np.flip(self.w.coefficients,axis=0)
        phi1 = np.polymul(w, np.polysub(np.polymul(u,u), np.polymul(v,v)))
        phi2 = np.polymul(2, np.polymul(w, np.polymul(u, v)))
        phi3 = np.polymul(1j, np.polymul(w, np.polyadd(np.polymul(u,u), np.polymul(v,v))))
        return phi1, phi2, phi3

    # 복소함수(계수로 주어짐)를 복소평면상의 0에서 해당 복소수까지 선적분(부정적분)하여 얻은 복소함수의 계수를 출력하는 함수.
    @staticmethod
    def poly_complex_integration(poly):
        n = len(list(poly))
        aux = np.array([1/n for n in np.arange(start=n,stop=0,step=-1)])
        out = aux*np.array(poly)
        return np.append(out,0)

    # Enneper-Weierstrass parameterization을 통해 얻은 함수로부터 평면을 나타내는 복소함수의 계수를 출력하는 함수.
    def get_surface_complex_poly_coef(self, mode='raw'): # 높은차수부터 낮은차수 순서로 출력함. 논문의 결과와 비교해서 검증완료.
        phi1_poly, phi2_poly, phi3_poly = self.get_holomorphic_coefs() 
        rx_poly = MSPN.poly_complex_integration(phi1_poly)
        ry_poly = MSPN.poly_complex_integration(phi2_poly)
        rz_poly = MSPN.poly_complex_integration(phi3_poly)
        if mode == "raw":
            return rx_poly, ry_poly, rz_poly
        elif mode == "numpy":
            r_coef_raw = (rx_poly, ry_poly, rz_poly)
            mlen = max([len(x) for x in r_coef_raw])
            r_coef = []
            for i in range(3):
                num_zeros = mlen - len(r_coef_raw[i])
                r_coef.append(np.concatenate([np.zeros(num_zeros), r_coef_raw[i]],axis=0))
            r_coef = np.array(r_coef)
            return r_coef
        else:
            raise ValueError(f'invalid mode argument: {mode}')

    # 평면을 나타내는 복소함수의 계수들로부터 평면의 함수를 출력하는 함수.
    def get_surface_handle(self, form='complex'):
        rx_poly, ry_poly, rz_poly = self.get_surface_complex_poly_coef()
        if form == 'real':
            def r(x, y):
                z = x + y*1j
                rx = np.real(np.polyval(rx_poly, z)) + self.r0[0]
                ry = np.real(np.polyval(ry_poly, z)) + self.r0[1]
                rz = np.real(np.polyval(rz_poly, z)) + self.r0[2]
                return rx, ry, rz
        elif form == "complex":
            def r(z):
                rx = np.real(np.polyval(rx_poly, z)) + self.r0[0]
                ry = np.real(np.polyval(ry_poly, z)) + self.r0[1]
                rz = np.real(np.polyval(rz_poly, z)) + self.r0[2]
                return rx, ry, rz
        else:
            raise ValueError('Invalid form argument')
        return r

    def symbolic(self): # Developing
        pass


class PlanarPHCurve():
    def __init__(self, HermiteInterpolation):
        self.H = HermiteInterpolation

    def __call__(self, tau, form="complex"):
        func = self.get_curve_handle(form=form)
        return func(tau)
        
    def get_curve_handle(self, form="complex"):
        if form == "complex":
            def curve(tau):
                q_real, q_imag = self.H.eval(tau)
                return q_real + q_imag*1j
        elif form == "real":
            def curve(tau):
                return self.H.eval(tau)
        return curve

    def get_curve_coef(self):
        return np.flip(PlanarPHCurve.BernsteinPoly_to_power_coef(self.H.cps),axis=1) # 높은차수부터 낮은차수 순서로 출력한다.

    @staticmethod
    def Bernsteinbasis_to_power_coef(n, k):
        coef = []
        for i in range(0,n+1):
            if i < k:
                coef.append(0)
            else:
                coef.append((-1)**(i-k)*math.comb(n,i)*math.comb(i,k))
        return coef

    @staticmethod
    def BernsteinPoly_to_power_coef(CPs): # 낮은차수의 계수부터 높은차수의 계수 순서로 출력한다.
        shape = np.shape(CPs) # [n-1, dim]
        n   = shape[0] - 1
        dim = shape[1]
        coef = None
        for k in range(n+1):
            power_coef = np.array(PlanarPHCurve.Bernsteinbasis_to_power_coef(n, k))
            if coef is None:
                coef = np.array([power_coef*cp for cp in CPs[k,:]])
            else:
                coef = coef + np.array([power_coef*cp for cp in CPs[k,:]])
        return coef


class SpacePHCurve():
    def __init__(self, MSPN, planarPHcurve):
        self.r = MSPN
        self.q = planarPHcurve

    def __call__(self, tau):
        func = self.get_curve_handle()
        return func(tau)

    def get_curve_handle(self):
        planar_curve_fun = self.q.get_curve_handle(form="complex")  
        surface_fun      = self.r.get_surface_handle(form="complex")
        def curve(tau):
            q_eval = planar_curve_fun(tau)
            r_eval = surface_fun(q_eval)
            return r_eval
        return curve

    def get_curve_poly_coef(self):
        q_coef = np.flip(self.q.get_curve_coef(), axis=1)                               # 높은차수부터 낮은차수 순서로 출력되므로 뒤집어 줌
        r_coef = np.flip(self.r.get_surface_complex_poly_coef(mode="numpy"), axis=1)    # 높은차수부터 낮은차수 순서로 출력되므로 뒤집어 줌
        return SpacePHCurve.compute_coef_of_composite(q_coef, r_coef)


    # 합성함수의 계수를 도출하는 코드 작성. 계수는 낮은 차수부터 높은 차수 순서로 쓴다.
    @staticmethod
    def compute_coef_of_composite(coef_q, coef_r): # 잘못됨. 버그 수정 필요.
        shape_q = np.shape(coef_q) # [2 x n+1]
        shape_r = np.shape(coef_r) # [3 x m+1]
        coef    = np.zeros([3,shape_q[1]+shape_r[1]], dtype=complex)
        for j in range(shape_r[1]):
            for k in range(shape_q[1]):
                power = j+k
                coef[0,power] += coef_r[0,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
                coef[1,power] += coef_r[1,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
                coef[2,power] += coef_r[2,j]*(coef_q[0,k] + coef_q[1,k]*1j)**j
        return coef


class SurfaceOptimizer():
    def __init__(self, wc, cc, w = [0.4, 0.3, 0.3]):
        self.WaypointConstraintsObj   = wc
        self.SpaceCurveConstraintsObj = cc
        self.weights                  = w

        self.WaypointConstraintsObj.weights = self.weights[0:2]

    def cost(self, X):
        return self.WaypointConstraintsObj(X) + self.SpaceCurveConstraintsObj(X)*self.weights[-1]

    def optimize(self, X0):
        res = scipy.optimize.minimize(self.cost, x0=X0)
        print(res['success'])
        Xopt = res['x']
        opt_MSPN = SurfaceOptimizer.MSPN_from_X(Xopt)
        return res, opt_MSPN

    def check_constraints(self, X):
        self.WaypointConstraintsObj(X, verbose=True)
        self.SpaceCurveConstraintsObj(X, verbose=True)
        return None

    @staticmethod
    def X2Params(X):
        r0            = X[:3]                                # constant term of the space curve
        preimage_coef = np.reshape(X[3:], newshape=[4,-1])   # poly. coefficients of complex functions, Phi.
        u_coef_real = preimage_coef[0,:]
        u_coef_imag = preimage_coef[1,:]
        v_coef_real = preimage_coef[2,:]
        v_coef_imag = preimage_coef[3,:]
        u_coef = u_coef_real + u_coef_imag*1j
        v_coef = v_coef_real + v_coef_imag*1j
        return r0, u_coef, v_coef

    @staticmethod
    def initialize_params(n):
        if n < 1:
            raise ValueError('n should be larger than 0.')
        nparam = 4*n + 7
        X0 = np.random.randn(nparam)
        return X0

    @staticmethod
    def MSPN_from_X(X):
        r0, u_coef, v_coef = SurfaceOptimizer.X2Params(X)
        u = PreimageComplexPolynomial(u_coef)
        v = PreimageComplexPolynomial(v_coef)
        w = PreimageComplexPolynomial([1+0j])
        r = MSPN(u, v, w, r0)
        return r

    def eval_cost(self):
        pass

    def compute_cost(self, X):
        pass


class WaypointConstraints():
    def __init__(self, Ps, Ts, qs, w=[0.5, 0.5]):
        self.space_waypoints_pos = Ps # shape = nwps+1 x 3
        self.space_waypoints_tan = Ts # shape = nwps+1 x 3
        self.planar_PHcurves     = qs # shape = nwps
        self.weights             = np.array(w)

    def __call__(self, X, verbose=False):
        pos_error_norm = self.eval_pos_error(X)
        tan_error_norm = self.eval_tan_error(X)
        if verbose:
            print(f"> Position error norm: \t{pos_error_norm}")
            print(f"> Tangent error norm: \t{tan_error_norm}")
        return self.weights[0]*pos_error_norm + self.weights[1]*tan_error_norm

    def eval_pos_error(self, X):
        nwps          = np.shape(self.space_waypoints_pos)[0]-1
        r             = SurfaceOptimizer.MSPN_from_X(X)
        pos_err_array = []
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            pos_error = np.array(p(0)) - np.array(self.space_waypoints_pos[i])
            pos_err_array.append(pos_error)
        error_norm_array = np.linalg.norm(pos_err_array,axis=1)
        return np.sum(error_norm_array)

    def eval_tan_error(self, X):
        nwps = np.shape(self.space_waypoints_pos)[0]-1
        r    = SurfaceOptimizer.MSPN_from_X(X)
        G    = []
        E    = []
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            grad = np.array(p(0.001)) - np.array(p(0))
            G.append(grad)
        G_u, T_u = [v / np.linalg.norm(v) for v in [G, self.space_waypoints_tan[:-1]]]
        for i in range(nwps):
            ang_btw  = np.arccos(np.clip(np.dot(G_u[i], T_u[i]), -1.0, 1.0))
            E.append(ang_btw)   # array of directional error (rad)
        return np.sum(E)


class SpaceCurveConstraints():
    def __init__(self, Ld, qs):
        self.length_desired = Ld
        self.planar_PHcurves= qs # shape = nwps

    def __call__(self, X, verbose=False):
        L  = self.eval_length(X)
        Le = np.abs(L-self.length_desired)
        if verbose:
            print(f"> Length error: \t{Le}   (Desired / Actual = {self.length_desired}, {L})")
        return Le

    def eval_length(self, X):
        nwps = np.shape(self.planar_PHcurves)[0]
        r    = SurfaceOptimizer.MSPN_from_X(X)
        Ls   = []
        tau  = np.linspace(0,1,100)
        for i in range(nwps):
            q = self.planar_PHcurves[i]
            p = SpacePHCurve(r, q)
            p_eval = p(tau)
            Ls.append(SpaceCurveConstraints.compute_length(p_eval))
        L = np.sum(Ls)
        return L
        
    # 3차원 곡선의 길이 수치적분하는 함수
    @staticmethod
    def compute_length(P):
        P = np.reshape(P, newshape=[3,-1])
        T = np.diff(P, axis=1)
        T = np.concatenate([T,np.reshape(T[:,0],newshape=[3,1])],axis=1)
        N = np.linalg.norm(T,axis=0)
        L = sum(N)
        return L
    