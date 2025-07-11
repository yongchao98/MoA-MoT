import numpy as np
from scipy.integrate import dblquad

def solve_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by computing the
    linking number of the preimages of two points on the sphere.
    """
    # 1. Define parameters based on the problem description.
    # The preimage of the South Pole (n_z = -1) is a circle C_S in the z=0 plane.
    # Its radius is sqrt(0.5).
    RHO_S = np.sqrt(0.5)
    
    # The preimage of the equatorial point (nx=1, ny=0, nz=0) is a curve C_E.
    # The condition nz=0 gives G = PI/2, which in turn defines a constant value
    # for r2. This value acts as a radius in a transformed coordinate system.
    # r2 = ln(2)/10
    R = np.log(2) / 10.0
    PI = np.pi

    print("Step 1: The Hopf charge is calculated as the linking number of two curves.")
    print(f"   - Curve C_S (preimage of South Pole) is a circle of radius {RHO_S:.4f} in the z=0 plane.")
    print(f"   - Curve C_E (preimage of an equatorial point) is a loop in the y=0 half-plane.")
    print("\nStep 2: The linking number is computed using the Gauss double integral formula.")
    print("   H = (1 / 4*pi) * Integral( (r1-r2) . (dr1 x dr2) / |r1-r2|^3 )")

    # 2. Define the integrand for the Gauss linking formula.
    # The function takes integration variables s (for C_E) and t (for C_S).
    def integrand(s, t, R_param, RHO_S_param):
        # Parametrization of C_S(t)
        rs_vec = np.array([RHO_S_param * np.cos(t), RHO_S_param * np.sin(t), 0])
        # Derivative of C_S(t) w.r.t. t
        drs_dt_vec = np.array([-RHO_S_param * np.sin(t), RHO_S_param * np.cos(t), 0])
    
        # Parametrization of C_E(s)
        x_e_sq = 0.5 + R_param * np.cos(s)
        x_e = np.sqrt(x_e_sq)
        re_vec = np.array([x_e, 0, R_param * np.sin(s)])
    
        # Derivative of C_E(s) w.r.t. s
        dx_e_ds = -R_param * np.sin(s) / (2.0 * x_e)
        dre_ds_vec = np.array([dx_e_ds, 0, R_param * np.cos(s)])
    
        # Numerator and denominator of the Gauss formula
        r_diff = rs_vec - re_vec
        cross_prod = np.cross(drs_dt_vec, dre_ds_vec)
        dot_prod = np.dot(r_diff, cross_prod)
        r_dist_cubed = np.sum(r_diff**2)**1.5
    
        # The curves do not intersect, so r_dist_cubed is never zero.
        return dot_prod / r_dist_cubed

    # 3. Perform the numerical integration.
    # The integral is over t from 0 to 2*pi (for C_S) and s from 0 to 2*pi (for C_E).
    # The additional parameters R and RHO_S are passed using the `args` tuple.
    integral_val, error = dblquad(integrand, 0, 2 * PI, lambda t: 0, lambda t: 2 * PI, args=(R, RHO_S))
    
    print("\nStep 3: Evaluate the integral and calculate the final charge.")
    print("The final equation for the Hopf charge H is:")
    print(f"   H = Integral Value / (4 * pi)")
    print(f"   H = {integral_val:.6f} / (4 * {PI:.6f})")

    # 4. Calculate the final Hopf charge.
    hopf_charge = integral_val / (4.0 * PI)
    print(f"   H ≈ {hopf_charge:.6f}")

    # The Hopf charge must be an integer. We round the numerical result.
    final_hopf_charge = round(hopf_charge)
    print(f"\nRounding to the nearest integer, the Hopf charge is: {final_hopf_charge}")

solve_hopf_charge()