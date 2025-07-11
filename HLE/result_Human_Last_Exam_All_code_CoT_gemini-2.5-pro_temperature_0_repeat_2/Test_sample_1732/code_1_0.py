import math

def solve_and_check_d():
    """
    This function solves the problem by analyzing the two possible cases
    derived from the problem's conditions.
    """
    print("Let the arithmetic sequence be a_n = a_1 + (n-1)d.")
    print(f"Let the second sequence be b_n = (n^2 + n) / a_n.")
    print("Since {b_n} is an arithmetic sequence, we have 2*b_2 = b_1 + b_3.")
    print("Substituting the definitions leads to the equation: a_1^2 - 3*a_1*d + 2*d^2 = 0.")
    print("Factoring this gives (a_1 - d)(a_1 - 2d) = 0. So, we have two cases: a_1 = d or a_1 = 2d.\n")

    # --- Case 1: a_1 = d ---
    print("--- Case 1: a_1 = d ---")
    print("a_n = d + (n-1)d = nd")
    print("b_n = n(n+1) / (nd) = (n+1)/d")
    print("S_99 is the sum of a_n from n=1 to 99: S_99 = d * sum(n) = d * (99*100/2) = 4950d")
    print("T_99 is the sum of b_n from n=1 to 99: T_99 = (1/d) * sum(n+1) = (1/d) * (99*100/2 + 99) = 5049/d")
    print("The condition S_99 - T_99 = 99 becomes: 4950d - 5049/d = 99")
    print("Multiplying by d and rearranging gives a quadratic equation for d.")
    print("Dividing by 99, we get: 50d^2 - d - 51 = 0")
    
    a1, b1, c1 = 50, -1, -51
    print(f"Solving the equation: {a1}d^2 + ({b1})d + ({c1}) = 0")
    
    delta1 = b1**2 - 4*a1*c1
    d1_1 = (-b1 + math.sqrt(delta1)) / (2*a1)
    d1_2 = (-b1 - math.sqrt(delta1)) / (2*a1)
    
    print(f"The solutions for d are {d1_1} and {d1_2}.")
    
    valid_solutions = []
    if d1_1 > 1:
        print(f"Solution d = {d1_1} is valid since it is > 1.")
        valid_solutions.append(d1_1)
    else:
        print(f"Solution d = {d1_1} is not valid.")
        
    if d1_2 > 1:
        print(f"Solution d = {d1_2} is valid since it is > 1.")
        valid_solutions.append(d1_2)
    else:
        print(f"Solution d = {d1_2} is not valid.")
    print("-" * 20)

    # --- Case 2: a_1 = 2d ---
    print("\n--- Case 2: a_1 = 2d ---")
    print("a_n = 2d + (n-1)d = (n+1)d")
    print("b_n = n(n+1) / ((n+1)d) = n/d")
    print("S_99 is the sum of a_n from n=1 to 99: S_99 = d * sum(n+1) = d * 5049 = 5049d")
    print("T_99 is the sum of b_n from n=1 to 99: T_99 = (1/d) * sum(n) = (1/d) * 4950 = 4950/d")
    print("The condition S_99 - T_99 = 99 becomes: 5049d - 4950/d = 99")
    print("Multiplying by d and rearranging gives a quadratic equation for d.")
    print("Dividing by 99, we get: 51d^2 - d - 50 = 0")

    a2, b2, c2 = 51, -1, -50
    print(f"Solving the equation: {a2}d^2 + ({b2})d + ({c2}) = 0")

    delta2 = b2**2 - 4*a2*c2
    d2_1 = (-b2 + math.sqrt(delta2)) / (2*a2)
    d2_2 = (-b2 - math.sqrt(delta2)) / (2*a2)

    print(f"The solutions for d are {d2_1} and {d2_2}.")

    if d2_1 > 1:
        print(f"Solution d = {d2_1} is valid since it is > 1.")
        valid_solutions.append(d2_1)
    else:
        print(f"Solution d = {d2_1} is not valid.")
        
    if d2_2 > 1:
        print(f"Solution d = {d2_2} is valid since it is > 1.")
        valid_solutions.append(d2_2)
    else:
        print(f"Solution d = {d2_2} is not valid.")
    print("-" * 20)

    # --- Conclusion ---
    if len(valid_solutions) == 1:
        print(f"\nThe only valid solution is d = {valid_solutions[0]}")
    else:
        print(f"\nThe valid solutions are: {valid_solutions}")

solve_and_check_d()