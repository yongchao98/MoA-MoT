import numpy as np

def solve_particle_emitter_problem():
    """
    Calculates the minimum ratio of the cube of the surface area to the square of the volume
    of the region traversed by particles from an emitter.

    The problem is solved analytically by:
    1. Defining the volume (V) and surface area (S) of the traversed region, which is a paraboloid segment.
    2. Expressing V and S in terms of a dimensionless parameter E = v^2/(2gh).
    3. The ratio Y = S^3/V^2 becomes a function Y(E).
    4. Differentiating Y(E) with respect to E and setting it to 0 leads to a cubic equation for E:
       73*E^3 + 51*E^2 - 21*E + 1 = 0.
    5. This equation has one physically meaningful positive root, which corresponds to the minimum of the ratio.
       This root is E_0 = (11 - 4*sqrt(3))/73.
    6. This script calculates E_0 and substitutes it back into the expression for Y to find the minimum value.
    """

    # Step 1: Find the value of E that minimizes the ratio.
    # The derivative of the ratio function being zero simplifies to the cubic equation:
    # 73*E^3 + 51*E^2 - 21*E + 1 = 0.
    # We solve the quadratic equation 73*E^2 - 22*E + 1 = 0, which results from factoring out (E+1).
    a, b, c = 73, -22, 1
    discriminant = b**2 - 4 * a * c
    
    # There are two potential roots for E.
    E_sol1 = ( -b + np.sqrt(discriminant) ) / (2 * a)
    E_sol2 = ( -b - np.sqrt(discriminant) ) / (2 * a)
    
    # Further analysis shows that only one of these roots is valid.
    # The condition is that E^(1/2)*(3-5E) and (1+2E)^(1/2)*(1-7E) must have the same sign.
    # E_sol1 (~0.246) is an extraneous solution introduced by squaring.
    # E_sol2 (~0.0557) is the correct root.
    E_min = E_sol2
    
    # Step 2: Define the components of the ratio in terms of E.
    # The ratio Y = S^3 / V^2 can be written as Y = (16*pi/27) * f(E), where
    # f(E) = s_norm(E)^3 / v_norm(E)^2.
    
    # Normalized volume squared component, v_norm(E) = E*(1+E)^2
    v_norm_sq = (E_min * (1 + E_min)**2)**2
    
    # Normalized surface area component, s_norm(E)
    # s_norm = 3*E + E^2 + 2*E^(1/2)*(1+2*E)^(3/2)
    s_norm = (3 * E_min + E_min**2 + 
              2 * np.sqrt(E_min) * (1 + 2 * E_min)**1.5)
              
    s_norm_cubed = s_norm**3
    
    # The dimensionless part of the ratio
    f_min = s_norm_cubed / v_norm_sq
    
    # The final ratio including the constant factor
    min_ratio = (16 * np.pi / 27) * f_min

    # Print the final numerical result
    print("The minimum ratio of the cube of the surface area to the square of the volume is:")
    print(min_ratio)

solve_particle_emitter_problem()