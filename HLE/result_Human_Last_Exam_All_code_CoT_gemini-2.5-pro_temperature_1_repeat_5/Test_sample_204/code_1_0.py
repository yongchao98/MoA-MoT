import numpy as np
from scipy.integrate import quad

def calculate_hopf_charge():
    """
    Calculates the Hopf charge by numerically computing the Gauss linking integral
    for the preimages of two points on the sphere.
    """
    PI = np.pi
    # Radius for the second curve C2, derived from G = PI/2
    R = np.log(2) / 10.0

    # Parameterization for curve C1 (preimage of the South Pole)
    # C1: r1(t) = (sqrt(0.5)*cos(t), sqrt(0.5)*sin(t), 0)
    def r1(t):
        return np.array([np.sqrt(0.5) * np.cos(t), np.sqrt(0.5) * np.sin(t), 0.0])

    def dr1_dt(t):
        return np.array([-np.sqrt(0.5) * np.sin(t), np.sqrt(0.5) * np.cos(t), 0.0])

    # Parameterization for curve C2 (preimage of an equatorial point)
    # C2: r2(s) = (0, sqrt(0.5 + R*cos(s)), R*sin(s))
    def r2(s):
        return np.array([0.0, np.sqrt(0.5 + R * np.cos(s)), R * np.sin(s)])

    def dr2_ds(s):
        # Derivative of y2(s) w.r.t. s
        dy2_ds = 0.5 * (1.0 / np.sqrt(0.5 + R * np.cos(s))) * (-R * np.sin(s))
        # Derivative of z2(s) w.r.t. s
        dz2_ds = R * np.cos(s)
        return np.array([0.0, dy2_ds, dz2_ds])

    # Integrand of the Gauss linking integral
    # The inner integral is over s (from 0 to 2*pi) for a fixed t
    def inner_integrand(s, t):
        vec_r1 = r1(t)
        vec_dr1 = dr1_dt(t)
        vec_r2 = r2(s)
        vec_dr2 = dr2_ds(s)

        r1_minus_r2 = vec_r1 - vec_r2
        r1_minus_r2_norm = np.linalg.norm(r1_minus_r2)

        # Avoid division by zero if curves touch (they don't in this case)
        if r1_minus_r2_norm < 1e-9:
            return 0.0

        # Numerator: (r1-r2) . (dr1 x dr2)
        cross_product = np.cross(vec_dr1, vec_dr2)
        numerator = np.dot(r1_minus_r2, cross_product)
        denominator = r1_minus_r2_norm**3
        
        return numerator / denominator

    # Function for the outer integral, which integrates the result of the inner one
    def outer_integrand(t):
        # Integrate inner_integrand over s from 0 to 2*pi
        integral_s, _ = quad(inner_integrand, 0, 2 * PI, args=(t,))
        return integral_s

    # Perform the outer integral over t from 0 to 2*pi
    linking_integral, _ = quad(outer_integrand, 0, 2 * PI)

    # Hopf charge is the linking number
    hopf_charge = linking_integral / (4 * PI)

    print("The final equation is: Hopf Charge = (1 / (4 * PI)) * I")
    print(f"where the numerical value of the Gauss linking integral I is: {linking_integral:.4f}")
    print(f"Hopf Charge = (1 / (4 * {PI:.4f})) * {linking_integral:.4f}")
    print(f"The calculated Hopf charge is: {hopf_charge:.4f}")
    
    # The Hopf charge is an integer, so we round the result.
    final_answer = int(round(hopf_charge))
    print(f"The integer Hopf charge is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    charge = calculate_hopf_charge()
    # The final answer is submitted in the specified format.
    # print(f'<<<{charge}>>>')