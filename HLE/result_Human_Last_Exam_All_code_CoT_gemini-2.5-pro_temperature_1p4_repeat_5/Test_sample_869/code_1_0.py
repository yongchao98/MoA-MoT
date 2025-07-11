import math

def solve_probability():
    """
    Calculates the probability based on the plan described.
    """
    f = math.factorial

    # Define constants from the problem
    N = 25
    T = 5
    items_per_type = 5
    individuals = 5
    items_per_individual = 5

    # Factorial terms
    f_N = f(N)
    f_items_per_type = f(items_per_type)
    f_individuals = f(individuals)
    f_2 = f(2)

    # 1. Calculate the total number of ways (S)
    # S = 25! / (5!)^5
    S_numerator = f_N
    S_denominator = f_items_per_type**T
    S = S_numerator / S_denominator

    print(f"Total number of distributions (S):")
    print(f"S = {N}! / ({items_per_type}!)**{T}")
    print(f"S = {S_numerator} / {S_denominator}")
    print(f"S = {S}\n")

    # 2. Calculate the number of favorable ways (F)
    # Based on the assumption of a symmetric distribution matrix C with C_ii = 2.
    # The number of such configurations is 6.
    num_configs = 6
    # The number of sequences for one such configuration N(C) = (5!)^5 / (2!)^5
    # The number of favorable sequences F = 5! * num_configs * N(C)
    F_numerator = f_individuals * num_configs * (f_items_per_type**individuals)
    F_denominator = f_2**individuals
    F = F_numerator / F_denominator
    
    print(f"Number of favorable distributions (F):")
    print(f"F = {individuals}! * {num_configs} * ({items_per_type}!)**{individuals} / ({2}!)**{individuals}")
    print(f"F = {f_individuals} * {num_configs} * {f_items_per_type**individuals} / {f_2**individuals}")
    print(f"F = {F}\n")

    # 3. Calculate the probability P = F / S
    # P = [5! * 6 * (5!)^5 / (2!)^5] / [25! / (5!)^5]
    # P = (5! * 6 * (5!)^10) / (25! * (2!)^5)
    P = F / S

    print(f"The probability P = F / S:")
    print(f"P = {F} / {S}")
    print(f"P = {P}\n")
    
    # Print the final equation with numbers as requested
    final_eq_num = f_individuals * num_configs * (f_items_per_type)**(T*2)
    final_eq_den = f_N * (f_2)**individuals
    
    print("Final Equation with Values:")
    print(f"P = ({f_individuals} * {num_configs} * ({f_items_per_type})**(10)) / ({f_N} * ({f_2})**{individuals})")
    print(f"P = {final_eq_num} / {final_eq_den}")
    
solve_probability()
<<<0.0009130741346393529>>>