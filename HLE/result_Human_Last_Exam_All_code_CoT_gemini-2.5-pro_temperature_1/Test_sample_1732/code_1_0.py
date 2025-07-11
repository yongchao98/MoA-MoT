import math

def solve_problem():
    """
    Solves for the common difference 'd' based on the problem conditions.
    """
    n = 99
    
    print("This script solves for the common difference 'd' of an arithmetic sequence.")
    print("Based on the derivation, there are two possible cases for the first term a_1: a_1 = d or a_1 = 2d.")
    print("\n" + "="*25)
    print("Evaluating Case 1: a_1 = d")
    print("="*25)
    
    # Case 1: a_1 = d
    # S_n = d * n(n+1)/2
    s_n_coeff_1 = n * (n + 1) / 2
    # T_n = (1/d) * (n(n+1)/2 + n)
    t_n_num_1 = s_n_coeff_1 + n
    
    print(f"For n = {n}:")
    print(f"S_{n} = {s_n_coeff_1}d")
    print(f"T_{n} = {t_n_num_1}/d")
    
    # Equation: s_n_coeff_1 * d - t_n_num_1 / d = n
    # s_n_coeff_1 * d^2 - n * d - t_n_num_1 = 0
    # Divide by n: (s_n_coeff_1/n) * d^2 - d - (t_n_num_1/n) = 0
    a1 = s_n_coeff_1 / n
    b1 = -1
    c1 = -t_n_num_1 / n
    
    print(f"\nThe condition S_{n} - T_{n} = {n} leads to the quadratic equation:")
    print(f"{int(a1)}d^2 + ({int(b1)})d + ({int(c1)}) = 0")

    # Solve the quadratic equation for d
    delta1 = b1**2 - 4 * a1 * c1
    d1_sol1 = (-b1 + math.sqrt(delta1)) / (2 * a1)
    d1_sol2 = (-b1 - math.sqrt(delta1)) / (2 * a1)
    
    print(f"The solutions for d are {d1_sol1:.2f} and {d1_sol2:.2f}.")
    
    valid_d = None
    if d1_sol1 > 1:
        valid_d = d1_sol1
        print(f"The solution d = {d1_sol1:.2f} satisfies the condition d > 1.")
    elif d1_sol2 > 1:
        valid_d = d1_sol2
        print(f"The solution d = {d1_sol2:.2f} satisfies the condition d > 1.")
    else:
        print("No solution in this case satisfies d > 1.")

    print("\n" + "="*25)
    print("Evaluating Case 2: a_1 = 2d")
    print("="*25)
    
    # Case 2: a_1 = 2d
    # S_n = d * (n(n+1)/2 + n)
    s_n_coeff_2 = t_n_num_1
    # T_n = (1/d) * n(n+1)/2
    t_n_num_2 = s_n_coeff_1
    
    print(f"For n = {n}:")
    print(f"S_{n} = {s_n_coeff_2}d")
    print(f"T_{n} = {t_n_num_2}/d")
    
    # Equation is s_n_coeff_2 * d - t_n_num_2 / d = n
    # Simplified: (s_n_coeff_2/n) * d^2 - d - (t_n_num_2/n) = 0
    a2 = s_n_coeff_2 / n
    b2 = -1
    c2 = -t_n_num_2 / n
    
    print(f"\nThe condition S_{n} - T_{n} = {n} leads to the quadratic equation:")
    print(f"{int(a2)}d^2 + ({int(b2)})d + ({int(c2)}) = 0")
    
    # Solve the quadratic equation
    delta2 = b2**2 - 4 * a2 * c2
    d2_sol1 = (-b2 + math.sqrt(delta2)) / (2 * a2)
    d2_sol2 = (-b2 - math.sqrt(delta2)) / (2 * a2)
    
    print(f"The solutions for d are {d2_sol1:.2f} and {d2_sol2:.2f}.")
    print("Neither solution satisfies the condition d > 1.")
    
    print("\n" + "="*25)
    print("Conclusion")
    print("="*25)
    if valid_d is not None:
        print(f"The only valid solution for d is {valid_d:.2f}.")
    else:
        print("No valid solution was found.")

solve_problem()