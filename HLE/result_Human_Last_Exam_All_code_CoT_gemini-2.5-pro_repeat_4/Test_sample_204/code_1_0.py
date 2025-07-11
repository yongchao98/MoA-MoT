import numpy as np
from scipy.integrate import dblquad

def solve_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by computing the
    linking number of the preimages of two points on the sphere S^2.
    """
    pi = np.pi

    # The vector field topology implies the Hopf charge is the linking number
    # of the preimages of two points. We choose the south pole and an
    # equatorial point.

    # 1. Define the preimage curves.
    # Preimage of the South Pole (nz=-1) is a circle C_S in the xy-plane.
    # G = pi => r2 = 0 => (x*x+y*y-0.5)^2 + z*z = 0 => z=0, x*x+y*y=0.5
    r_S_radius = np.sqrt(0.5)

    # Preimage of an equatorial point (nx=1, ny=0, nz=0) is a loop C_E.
    # nz=0 => G=pi/2 => pi*exp(-10*r2)=pi/2 => r2=ln(2)/10
    # (x*x+y*y-0.5)^2 + z*z = (ln(2)/10)^2
    # f=0 => y=0, x>0. This defines a loop C_E in the xz-plane.
    r_E_torus_radius = np.log(2) / 10.0

    # 2. Parameterize the curves for the Gauss linking integral.
    # Curve C_S(t) = (r_S*cos(t), r_S*sin(t), 0)
    def r_S_vec(t):
        return np.array([r_S_radius * np.cos(t), r_S_radius * np.sin(t), 0.0])

    def dr_S_dt_vec(t):
        return np.array([-r_S_radius * np.sin(t), r_S_radius * np.cos(t), 0.0])

    # Curve C_E(s) in the xz-plane
    # Let x^2-0.5 = r_E_torus_radius*cos(s), z = r_E_torus_radius*sin(s)
    # x = sqrt(0.5 + r_E_torus_radius*cos(s))
    def r_E_vec(s):
        x_E = np.sqrt(0.5 + r_E_torus_radius * np.cos(s))
        z_E = r_E_torus_radius * np.sin(s)
        return np.array([x_E, 0.0, z_E])

    def dr_E_ds_vec(s):
        # Using the chain rule for dx/ds
        cos_s = np.cos(s)
        sin_s = np.sin(s)
        x_E = np.sqrt(0.5 + r_E_torus_radius * cos_s)
        dx_ds = (-r_E_torus_radius * sin_s) / (2.0 * x_E)
        dz_ds = r_E_torus_radius * cos_s
        return np.array([dx_ds, 0.0, dz_ds])

    # 3. Define the integrand for the Gauss linking formula.
    # Integrand = (r_E(s) - r_S(t)) . (dr_E(s) x dr_S(t)) / |r_E(s) - r_S(t)|^3
    # Note the argument order for dblquad: func(y, x), so we use (t, s)
    def gauss_integrand(t, s):
        rE = r_E_vec(s)
        rS = r_S_vec(t)
        drE = dr_E_ds_vec(s)
        drS = dr_S_dt_vec(t)

        r_diff = rE - rS
        r_dist = np.linalg.norm(r_diff)
        
        # The curves do not intersect, so dist should not be zero.
        if r_dist < 1e-9:
            return 0.0

        cross_prod = np.cross(drE, drS)
        numerator = np.dot(r_diff, cross_prod)
        denominator = r_dist**3
        
        return numerator / denominator

    # 4. Compute the double integral numerically.
    # The integral is over s from 0 to 2*pi and t from 0 to 2*pi.
    print("Calculating the Gauss linking integral numerically...")
    integral_val, integral_err = dblquad(
        gauss_integrand, 0, 2 * pi, lambda s: 0, lambda s: 2 * pi
    )

    # 5. Calculate the Hopf charge from the integral value.
    # Hopf Charge Q = (1 / 4π) * integral_val
    four_pi = 4 * pi
    hopf_charge = integral_val / four_pi
    
    print(f"\nThe value of the double integral is: {integral_val:.4f}")
    print("The Hopf charge is given by the formula: Q = Integral / (4 * PI)")
    print(f"Q = {integral_val:.4f} / (4 * {pi:.4f})")
    print(f"Q = {integral_val:.4f} / {four_pi:.4f}")
    print(f"Calculated Hopf Charge = {hopf_charge:.4f}")
    
    # The Hopf charge is a topological integer.
    final_answer = round(hopf_charge)
    print(f"\nRounding to the nearest integer, the Hopf charge is {final_answer}.")
    return final_answer

if __name__ == '__main__':
    solve_hopf_charge()