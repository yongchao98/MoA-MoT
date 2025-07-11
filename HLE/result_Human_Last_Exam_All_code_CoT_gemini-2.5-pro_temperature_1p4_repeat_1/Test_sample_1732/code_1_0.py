import math

def solve_problem():
    """
    Solves for the common difference 'd' based on the problem's conditions.
    The derivation leads to two possible quadratic equations for d. This function
    solves both and finds the valid solution.
    """
    print("Based on the problem statement, we have two possible cases.")
    print("-" * 50)
    
    # --- Case 1: a_n = dn ---
    print("Case 1: a_n = dn")
    print("The condition S_99 - T_99 = 99 leads to the quadratic equation for d:")
    a1, b1, c1 = 50, -1, -51
    print(f"{a1} * d^2 + ({b1}) * d + ({c1}) = 0")
    
    # Calculate the discriminant
    delta1 = b1**2 - 4 * a1 * c1
    
    valid_d_found = False
    final_d = None
    
    if delta1 >= 0:
        # Calculate the two roots
        d1_root1 = (-b1 + math.sqrt(delta1)) / (2 * a1)
        d1_root2 = (-b1 - math.sqrt(delta1)) / (2 * a1)
        
        print(f"The solutions for d are {d1_root1:.2f} and {d1_root2:.2f}.")
        
        # Check the condition d > 1
        print("Checking which solution satisfies d > 1:")
        if d1_root1 > 1:
            print(f"{d1_root1:.2f} is a valid solution.")
            valid_d_found = True
            final_d = d1_root1
        else:
            print(f"{d1_root1:.2f} is not greater than 1.")
            
        if d1_root2 > 1:
            print(f"{d1_root2:.2f} is a valid solution.")
            valid_d_found = True
            final_d = d1_root2
        else:
            print(f"{d1_root2:.2f} is not greater than 1.")
    else:
        print("No real solutions for d in this case.")

    print("-" * 50)
    
    # --- Case 2: a_n = d(n+1) ---
    print("Case 2: a_n = d(n+1)")
    print("The condition S_99 - T_99 = 99 leads to the quadratic equation for d:")
    a2, b2, c2 = 51, -1, -50
    print(f"{a2} * d^2 + ({b2}) * d + ({c2}) = 0")
    
    # Calculate the discriminant
    delta2 = b2**2 - 4 * a2 * c2
    
    if delta2 >= 0:
        # Calculate the two roots
        d2_root1 = (-b2 + math.sqrt(delta2)) / (2 * a2)
        d2_root2 = (-b2 - math.sqrt(delta2)) / (2 * a2)

        print(f"The solutions for d are {d2_root1:.2f} and {d2_root2:.2f}.")
        
        # Check the condition d > 1
        print("Checking which solution satisfies d > 1:")
        if d2_root1 > 1:
            print(f"{d2_root1:.2f} is a valid solution.")
            valid_d_found = True
            final_d = d2_root1
        else:
            print(f"{d2_root1:.2f} is not greater than 1.")
            
        if d2_root2 > 1:
            print(f"{d2_root2:.2f} is a valid solution.")
            valid_d_found = True
            final_d = d2_root2
        else:
            print(f"{d2_root2:.2f} is not greater than 1.")
    else:
        print("No real solutions for d in this case.")
        
    print("-" * 50)
    
    if valid_d_found:
        print(f"The only solution that satisfies all conditions is d = {final_d:.2f}")
    else:
        print("No solution found that satisfies all conditions.")

# Run the solver
solve_problem()