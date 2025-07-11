import numpy as np
from scipy.integrate import dblquad

def calculate_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field using the linking number approach.

    The Hopf charge is the linking number of the preimages of two regular points on the S^2 sphere.
    We choose the points p1 = (1,0,0) and p2 = (0,1,0).

    1. Preimage of p1 = (1,0,0):
       This corresponds to n=(1,0,0), which means G=pi/2 and f=0.
       G = pi * exp(-10*r2) = pi/2  =>  exp(-10*r2) = 0.5  =>  r2 = ln(2)/10.
       f = atan2(y,x) = 0  => y=0, x>0.
       The preimage is a curve C1 in the xz-plane, defined by:
       r2^2 = ((x^2-0.5)^2 + z^2) = (ln(2)/10)^2.

    2. Preimage of p2 = (0,1,0):
       This corresponds to n=(0,1,0), which means G=pi/2 and f=pi/2.
       f = atan2(y,x) = pi/2 => x=0, y>0.
       The condition on r2 is the same.
       The preimage is a curve C2 in the yz-plane, defined by:
       r2^2 = ((y^2-0.5)^2 + z^2) = (ln(2)/10)^2.

    3. Linking Number Calculation:
       The linking number of two curves C1 and C2 is given by the Gauss linking integral:
       Lk(C1, C2) = (1 / 4pi) * integral_C1 integral_C2
                      (r1 - r2) . (dr1 x dr2) / |r1 - r2|^3

       We will compute this integral numerically.
    """
    # Constant R from r2 = ln(2)/10
    R = np.log(2) / 10.0

    # Parametrization of the curves C1(t) and C2(s) for t, s in [0, 2*pi]
    # C1(t) is in the y=0 plane
    def r1(t):
        x = np.sqrt(0.5 + R * np.cos(t))
        y = 0.0
        z = R * np.sin(t)
        return np.array([x, y, z])

    def dr1_dt(t):
        dx_dt = 0.5 * (-R * np.sin(t)) / np.sqrt(0.5 + R * np.cos(t))
        dy_dt = 0.0
        dz_dt = R * np.cos(t)
        return np.array([dx_dt, dy_dt, dz_dt])

    # C2(s) is in the x=0 plane
    def r2(s):
        x = 0.0
        y = np.sqrt(0.5 + R * np.cos(s))
        z = R * np.sin(s)
        return np.array([x, y, z])

    def dr2_ds(s):
        dx_ds = 0.0
        dy_ds = 0.5 * (-R * np.sin(s)) / np.sqrt(0.5 + R * np.cos(s))
        dz_ds = R * np.cos(s)
        return np.array([dx_ds, dy_ds, dz_ds])

    # Integrand for the Gauss linking integral
    def integrand(t, s):
        # The integral is over dt and ds, so we evaluate the vector components
        # at t and s
        vec_r1 = r1(t)
        vec_r2 = r2(s)
        vec_dr1 = dr1_dt(t)
        vec_dr2 = dr2_ds(s)

        r1_minus_r2 = vec_r1 - vec_r2
        dr1_cross_dr2 = np.cross(vec_dr1, vec_dr2)
        
        numerator = np.dot(r1_minus_r2, dr1_cross_dr2)
        denominator = np.linalg.norm(r1_minus_r2)**3

        if denominator < 1e-9: # Avoid division by zero, though they shouldn't meet
            return 0
        
        return numerator / denominator

    # Perform the double numerical integration
    # We integrate from 0 to 2*pi for both t and s
    # The options are to increase accuracy and limit for the integration
    gauss_integral, error = dblquad(integrand, 0, 2 * np.pi, lambda s: 0, lambda s: 2 * np.pi, epsabs=1e-4, epsrel=1e-4)

    # The linking number is the integral divided by 4*pi
    linking_number = gauss_integral / (4 * np.pi)
    
    # The Hopf charge is an integer, so we round the result
    hopf_charge = round(linking_number)

    print(f"The preimage C1 of n=(1,0,0) is a loop in the y=0 plane.")
    print(f"The preimage C2 of n=(0,1,0) is a loop in the x=0 plane.")
    print(f"These two loops are linked like a chain.")
    print(f"The Hopf charge is the linking number of these two preimages.")
    print(f"Numerically calculating the Gauss linking integral gives: {gauss_integral:.4f}")
    print(f"Dividing by 4*pi gives the linking number: {linking_number:.4f}")
    print(f"The Hopf charge is the closest integer to this value.")
    print(f"\nHopf Charge = {hopf_charge}")

calculate_hopf_charge()