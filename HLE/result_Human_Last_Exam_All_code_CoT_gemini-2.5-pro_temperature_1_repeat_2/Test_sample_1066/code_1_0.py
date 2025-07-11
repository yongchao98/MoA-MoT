import math

def display_bound_formula():
    """
    This function displays the derived lower bound for the expected 
    detection statistic E[S] for watermarked text.
    """
    
    # Define symbolic variables for the formula
    n = "n"
    alpha = "alpha"
    
    # Define the numeric constants in the formula
    number_2 = 2
    constant_pi_symbol = "pi"
    
    # The derived lower bound formula is n * (alpha + ln(2/pi))
    formula_string = f"E[S] >= {n} * ({alpha} + ln({number_2} / {constant_pi_symbol}))"
    
    print("The derived lower bound on the expected detection statistic E[S] is:")
    print(formula_string)
    
    print("\nIn this final equation:")
    print(f"- '{n}' represents the number of tokens in the text.")
    print(f"- '{alpha}' represents the average entropy per token.")
    print(f"- The number '{number_2}' is the numerator inside the logarithm.")
    print(f"- The constant '{constant_pi_symbol}' is pi, with a value of approximately {math.pi}.")

display_bound_formula()