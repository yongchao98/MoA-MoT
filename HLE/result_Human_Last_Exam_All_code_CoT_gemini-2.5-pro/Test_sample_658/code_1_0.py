import math

def print_asymptotic_formula():
    """
    Prints the derived asymptotic formula for d_{B,delta}.
    """
    # The asymptotic value is a formula A(B, delta).
    # We use L = log(delta^-1) in the formula.
    # The derived formula that combines the two regimes is B + L / log(L/B + C).
    # We choose the constant C = e.
    B = "B"
    L = "L"
    log_func = "log"
    constant = str(round(math.e, 3)) # The number in the equation

    # Construct the formula string
    formula = f"{B} + {L} / {log_func}({L}/{B} + {constant})"
    
    print(formula)

print_asymptotic_formula()