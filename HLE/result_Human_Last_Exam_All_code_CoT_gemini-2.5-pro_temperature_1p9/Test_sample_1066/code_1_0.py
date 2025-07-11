import math

def print_lower_bound_formula():
    """
    This function prints the lower bound formula for the expected detection score E[S].
    The formula involves the number of tokens n, the average entropy per token alpha,
    and the mathematical constant pi.
    """
    n_str = "n"
    alpha_str = "\u03B1"  # Unicode for alpha
    pi_str = "\u03C0"    # Unicode for pi
    
    # The derived lower bound for the expected score E[S] is n * (alpha + 1 - pi^2 / 12)
    # We will print each part of the formula.
    
    term_1 = 1
    term_2 = math.pi**2
    term_3 = 12
    
    final_constant = term_1 - term_2 / term_3
    
    print("A lower bound on the expected detection score E[S] is given by the formula:")
    print(f"E[S] >= {n_str} * ({alpha_str} + {term_1} - {pi_str}^2 / {term_3})")
    print("\nWhere:")
    print(f"  n = number of tokens in the text")
    print(f"  {alpha_str} = average entropy per token as perceived by the language model")
    print(f"  {pi_str} = {math.pi}")
    print("\nCalculating the constants:")
    print(f"  {pi_str}^2 = {term_2}")
    print(f"  The constant part of the formula per token is {term_1} - {term_2} / {term_3} \u2248 {final_constant}")
    
print_lower_bound_formula()