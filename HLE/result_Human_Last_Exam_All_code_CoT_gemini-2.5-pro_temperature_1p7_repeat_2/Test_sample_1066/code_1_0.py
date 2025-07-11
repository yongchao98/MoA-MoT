import math

def display_lower_bound_formula():
    """
    Prints the mathematical formula for the lower bound of the expected 
    watermark score E[S], as requested by the user.
    The formula involves the number of tokens (n), the average entropy (alpha),
    and the constant pi.
    """
    
    # Define the symbols for the equation. The prompt asks to output each number.
    n_symbol = "n"
    alpha_symbol = "α"
    pi_symbol = "π"
    constant_6 = 6
    constant_2 = 2
    
    # The lower bound is n*alpha + n*ln(6/pi^2).
    # We will print this expression in a formatted way.
    
    print("A lower bound on the expected detection score E[S] for a watermarked text is given by the formula:")
    print(f"E[S] >= {n_symbol} * {alpha_symbol} + {n_symbol} * ln({constant_6} / {pi_symbol}^{constant_2})")

# Execute the function to display the answer.
display_lower_bound_formula()
