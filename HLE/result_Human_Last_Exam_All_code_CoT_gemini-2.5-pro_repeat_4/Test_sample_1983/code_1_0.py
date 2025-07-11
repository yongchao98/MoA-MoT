import math

def solve_dichotomy_problem():
    """
    Solves the mathematical problem based on the theory of discrete dichotomies.
    """
    # 1. Define the parameters from the problem statement.
    # We interpret k1 and k2 from their log10 values to handle large numbers.
    log10_k1 = 3000
    log10_k2 = 500
    
    # Other parameters
    lambda1 = 0.5
    lambda2 = 0.5
    h_norm = 1000

    # 2. State the formulas for the asymptotic norms.
    # lim_{n->+inf} ||x_n|| = k1 * ||h|| / (1 - lambda1)
    # lim_{n->-inf} ||x_n|| = k2 * ||h|| / (1 - lambda2)
    
    # 3. Calculate the two parts of the expression using logarithms to handle large numbers.
    # The expression is E = 100 * log10(1/3 * lim_plus_inf) + 10 * log10(1/3 * lim_minus_inf)

    # Part 1: 100 * log10( (1/3) * (k1 * h_norm) / (1 - lambda1) )
    # log10_part1 = log10(1/3) + log10(k1) + log10(h_norm) - log10(1-lambda1)
    log10_term1 = math.log10(1/3) + log10_k1 + math.log10(h_norm) - math.log10(1 - lambda1)
    part1 = 100 * log10_term1

    # Part 2: 10 * log10( (1/3) * (k2 * h_norm) / (1 - lambda2) )
    # log10_part2 = log10(1/3) + log10(k2) + log10(h_norm) - log10(1-lambda2)
    log10_term2 = math.log10(1/3) + log10_k2 + math.log10(h_norm) - math.log10(1 - lambda2)
    part2 = 10 * log10_term2
    
    # The final result is the sum of the two parts.
    result = part1 + part2

    # 4. Print the final equation with the numbers, as requested.
    # We calculate the values inside the logs for printing.
    lim_plus_inf_val = (10**log10_k1 * h_norm) / (1-lambda1)
    lim_minus_inf_val = (10**log10_k2 * h_norm) / (1-lambda2)
    
    # Use string formatting to represent large numbers in scientific notation.
    print("The final equation with the computed numbers is:")
    print(f"100 * log10( (1/3) * {lim_plus_inf_val:.1e} ) + 10 * log10( (1/3) * {lim_minus_inf_val:.1e} )")
    
    print("\nThe final numerical result is:")
    print(result)

# Execute the function to get the answer.
solve_dichotomy_problem()