import math

def print_lower_bound_formula():
    """
    This function prints the derived lower bound for E[S].
    The bound involves the number of tokens n, the average entropy alpha,
    and the constant pi.
    """
    
    # The derived lower bound for E[S] is n*alpha - n*ln(pi^2 / 6).
    # We will construct this formula as a string and print it.
    
    n_str = "n"
    alpha_str = "\u03B1"  # Unicode for alpha
    pi_str = "\u03C0"    # Unicode for pi
    
    # Using names for constants for clarity
    constant_name = "ln(" + pi_str + "^2 / 6)"
    
    # Constructing the full equation string
    # E[S] >= n*alpha - n*ln(pi^2/6)
    equation = f"E[S] >= {n_str}\u22C5{alpha_str} - {n_str}\u22C5{constant_name}"
    
    print("The derived lower bound is:")
    print(equation)
    
    # The prompt requests to "output each number in the final equation".
    # Let's calculate the numerical value of the constant factor ln(pi^2/6).
    
    pi_val = math.pi
    constant_val = math.log(pi_val**2 / 6)
    
    print("\nWhere:")
    print(f"n = number of tokens")
    print(f"\u03B1 = average token entropy")
    print(f"\u03C0 \u2248 {pi_val:.6f}")
    print(f"The constant term ln(\u03C0^2 / 6) \u2248 {constant_val:.6f}")

if __name__ == '__main__':
    print_lower_bound_formula()
