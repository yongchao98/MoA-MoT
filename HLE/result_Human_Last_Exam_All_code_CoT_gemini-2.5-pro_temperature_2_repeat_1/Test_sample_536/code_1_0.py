import math

def solve_for_angle():
    """
    This function solves for the limiting value of the inner product <b_p, z_p>.
    
    The core idea is based on the geometric interpretation of the problem. In the
    high-dimensional limit, the three unit vectors h_p, b_p, and z_p are assumed
    to be coplanar. This linear dependence can be expressed by stating that their
    Gram determinant is zero.
    
    Let:
    c_hb = lim <h_p, b_p> = 0.9375
    c_hz = lim <h_p, z_p> = 0.9
    c_bz = lim <b_p, z_p> (the value we want to find)
    
    The Gram determinant equation for these three unit vectors is:
    | 1    c_hb  c_hz |
    | c_hb  1    c_bz | = 0
    | c_hz c_bz   1   |
    
    Expanding this determinant leads to the following quadratic equation for c_bz:
    c_bz^2 - 2 * c_hb * c_hz * c_bz + (c_hb^2 + c_hz^2 - 1) = 0
    
    We will now substitute the given values and solve for c_bz.
    """
    
    c_hb = 0.9375
    c_hz = 0.9
    
    # Coefficients of the quadratic equation a*x^2 + b*x + c = 0 for x = c_bz
    a = 1.0
    b = -2 * c_hb * c_hz
    c = c_hb**2 + c_hz**2 - 1
    
    print(f"Based on the coplanarity assumption, we solve the quadratic equation:")
    print(f"c_bz^2 - 2*({c_hb})*({c_hz})*c_bz + (({c_hb})^2 + ({c_hz})^2 - 1) = 0")
    print(f"Simplified equation: {a:.4f}*c_bz^2 + {b:.4f}*c_bz + {c:.8f} = 0")
    print("-" * 30)
    
    # Calculate the discriminant
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("The discriminant is negative, no real solutions exist.")
        return

    # Solve for the two roots
    root1 = (-b + math.sqrt(discriminant)) / (2 * a)
    root2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    # In RMT models of this type, the leading sample eigenvector (h_p) typically lies
    # geometrically between the population eigenvector (b_p) and the principal
    # direction of the noise perturbation. This corresponds to the angle between
    # b_p and z_p being the difference of the other angles, which leads to choosing
    # the smaller positive root. We select the smaller root as the physically plausible one.
    solution = min(root1, root2)
    
    # We must pick one solution. Let the angles be a_xy = arccos(c_xy)
    # The two solutions correspond to a_hz = a_hb + a_bz and a_hz = |a_hb - a_bz|
    # The latter is the more physically plausible configuration.
    # The smaller root corresponds to this case.
    
    print(f"The two possible solutions for <b_p, z_p> are: {root1:.6f} and {root2:.6f}")
    print("Choosing the physically plausible solution.")
    print(f"The final calculated value for lim <b_p, z_p> is: {solution}")

solve_for_angle()
<<<0.692067>>>