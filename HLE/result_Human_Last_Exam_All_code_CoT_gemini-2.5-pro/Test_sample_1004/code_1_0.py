def solve_river_problem():
    """
    This script calculates the total downstream drift of a boat based on the problem description.
    It derives the final formula by breaking down the problem into steps.
    """
    print("This program calculates the final distance between the boat's starting and returning positions.")
    print("Let L be the river width, v_0 be the max flow velocity, and v be the boat's speed relative to the water.")
    print("-" * 70)

    # Step 1: Define the river velocity profile
    print("Step 1: Define the river's flow velocity v_flow(y).")
    print("The flow is proportional to the distance from the nearest shore.")
    print("v_flow(y) = (2 * v_0 / L) * y,          for 0 <= y <= L/2")
    print("v_flow(y) = (2 * v_0 / L) * (L - y),   for L/2 < y <= L")
    print("-" * 70)

    # Step 2: Formulate the drift integral
    print("Step 2: Formulate the total downstream drift.")
    print("The total drift 'x_total' is the integral of the flow velocity over the travel time.")
    print("dt = dy / v, so dx = v_flow(y) * (dy / v).")
    print("x_total = x_outward + x_return")
    print("-" * 70)

    # Step 3: Calculate drift for the outward trip (y from 0 to 3L/4)
    print("Step 3: Calculate the outward drift, x_outward, from y=0 to y=3L/4.")
    print("We must split the integral at y = L/2 because the flow function changes.")

    # The result of integrating (2*v_0/L)*y from 0 to L/2 is (v_0*L)/4
    coef1_num, coef1_den = 1, 4
    print(f"  - The drift from y=0 to y=L/2 is: ({coef1_num}/{coef1_den}) * (v_0 * L / v)")

    # The result of integrating (2*v_0/L)*(L-y) from L/2 to 3L/4 is (3*v_0*L)/16
    coef2_num, coef2_den = 3, 16
    print(f"  - The drift from y=L/2 to y=3L/4 is: ({coef2_num}/{coef2_den}) * (v_0 * L / v)")

    # Summing the coefficients for the outward trip: 1/4 + 3/16 = 4/16 + 3/16 = 7/16
    outward_num, outward_den = 7, 16
    print(f"Total outward drift, x_outward = (1/4 + 3/16) * (v_0 * L / v) = ({outward_num}/{outward_den}) * (v_0 * L / v)")
    print("-" * 70)

    # Step 4: Calculate drift for the return trip (y from 3L/4 to 0)
    print("Step 4: Calculate the return drift, x_return, from y=3L/4 to y=0.")
    print("The time spent at each distance 'y' from the bank is the same as the outward trip.")
    print("Therefore, the drift for the return trip is identical to the outward trip.")
    print(f"Total return drift, x_return = ({outward_num}/{outward_den}) * (v_0 * L / v)")
    print("-" * 70)

    # Step 5: Calculate total drift
    print("Step 5: Calculate the total drift, x_total = x_outward + x_return.")
    total_drift_num = 2 * outward_num
    total_drift_den = outward_den
    print(f"x_total = ({outward_num}/{outward_den})*(v_0*L/v) + ({outward_num}/{outward_den})*(v_0*L/v)")
    print(f"x_total = 2 * ({outward_num}/{outward_den}) * (v_0 * L / v) = ({total_drift_num}/{total_drift_den}) * (v_0*L/v)")
    
    final_num, final_den = 7, 8
    print("\nSimplifying the fraction gives the final equation:")
    print(f"Total Downstream Distance = ({final_num}/{final_den}) * (v_0 * L / v)")

    # The final coefficient is 7/8
    final_coefficient = final_num / final_den
    return final_coefficient

# Execute the function and print the final numerical coefficient
final_answer = solve_river_problem()
print(f"\nThe numerical coefficient of the expression is {final_answer}.")

<<<0.875>>>