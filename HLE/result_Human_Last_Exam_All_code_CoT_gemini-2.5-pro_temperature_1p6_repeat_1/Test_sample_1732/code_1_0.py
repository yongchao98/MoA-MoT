import math

def solve_sequence_problem():
    """
    Solves the math problem by following the derived steps,
    analyzing the two possible cases, and finding the correct value for d.
    """
    print("Let the arithmetic sequence {a_n} have first term a_1 and common difference d.")
    print("We are given a_n = a_1 + (n-1)d and b_n = n(n+1) / a_n.")
    print("Since {b_n} is an arithmetic progression, we have b_1 + b_3 = 2 * b_2.\n")

    print("Step 1: Express b_1, b_2, b_3 in terms of a_1 and d.")
    print("b_1 = 2 / a_1")
    print("b_2 = 6 / (a_1 + d)")
    print("b_3 = 12 / (a_1 + 2d)\n")

    print("Step 2: From b_1 + b_3 = 2*b_2, we derive the relation between a_1 and d.")
    print("Equation: 2/a_1 + 12/(a_1 + 2d) = 12/(a_1 + d)")
    print("This simplifies to the quadratic equation: a_1^2 - 3*a_1*d + 2*d^2 = 0")
    print("Factoring gives: (a_1 - d) * (a_1 - 2d) = 0")
    print("This leads to two cases: Case 1 (a_1 = d) and Case 2 (a_1 = 2d).\n")

    # Define n
    n = 99

    print(f"--- Analyzing Case 1: a_1 = d ---")
    print("If a_1 = d, then a_n = nd and b_n = (n+1)/d.")
    # Calculate coefficients for S_99 and T_99
    s99_coeff_c1 = n * (n + 1) // 2
    t99_num_c1 = n * (n + 1) // 2 + n
    print(f"S_99 = d * (99*100/2) = {s99_coeff_c1}d")
    print(f"T_99 = (1/d) * (99*100/2 + 99) = {t99_num_c1}/d")
    print(f"Using the condition S_99 - T_99 = 99:")
    print(f"The equation is: {s99_coeff_c1}d - {t99_num_c1}/d = {n}")
    print(f"This leads to the quadratic equation: {s99_coeff_c1}d^2 - {n}d - {t99_num_c1} = 0")
    # After dividing by 99
    a1, b1, c1 = 50, -1, -51
    print(f"Dividing by {n}, we get: {a1}d^2 - {b1}d - {c1} = 0")
    # Solve quadratic equation
    discriminant = b1**2 - 4 * a1 * c1
    d_sol1 = (-b1 + math.sqrt(discriminant)) / (2 * a1)
    d_sol2 = (-b1 - math.sqrt(discriminant)) / (2 * a1)
    print(f"The solutions are d = {d_sol1} and d = {d_sol2}.")
    print(f"Since d > 1, the valid solution for this case is d = {d_sol1}.\n")

    print(f"--- Analyzing Case 2: a_1 = 2d ---")
    print("If a_1 = 2d, then a_n = (n+1)d and b_n = n/d.")
    s99_coeff_c2 = n * (n + 1) // 2 + n
    t99_num_c2 = n * (n + 1) // 2
    print(f"S_99 = d * (99*100/2 + 99) = {s99_coeff_c2}d")
    print(f"T_99 = (1/d) * (99*100/2) = {t99_num_c2}/d")
    print(f"Using the condition S_99 - T_99 = 99:")
    print(f"The equation is: {s99_coeff_c2}d - {t99_num_c2}/d = {n}")
    print(f"This leads to the quadratic equation: {s99_coeff_c2}d^2 - {n}d - {t99_num_c2} = 0")
    # After dividing by 99
    a2, b2, c2 = 51, -1, -50
    print(f"Dividing by {n}, we get: {a2}d^2 - {b2}d - {c2} = 0")
    # Solve quadratic equation (by factoring: (51d+50)(d-1)=0)
    d_sol3 = 1.0
    d_sol4 = -50/51
    print(f"The solutions are d = {d_sol3} and d = {d_sol4:.4f}.")
    print("Neither solution satisfies the condition d > 1.\n")

    print("--- Conclusion ---")
    print("The only case that yields a valid result is Case 1.")
    print(f"The final value of d is {d_sol1}.")

solve_sequence_problem()