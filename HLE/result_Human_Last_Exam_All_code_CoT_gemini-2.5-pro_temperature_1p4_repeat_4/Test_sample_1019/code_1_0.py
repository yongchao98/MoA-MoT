import math

def calculate_hamiltonicity_threshold(n, eta):
    """
    Calculates the d-threshold for Hamiltonicity for H_n U G(n, p).

    Args:
        n (int): The number of vertices in the graph.
        eta (float): The parameter defining the minimum degree d = n/2 - eta.
    """

    # Check if inputs are valid for log calculations
    if n <= 1 or math.log(n) <= 0:
        print("Error: 'n' must be large enough for log(log(n)) to be positive (n > e).")
        return

    # Check if eta is in the specified range
    if not (0.5 <= eta <= n / 64):
        print(f"Warning: eta={eta} is outside the specified range [0.5, {n/64:.2f}] for n={n}.")


    # The d-threshold formula
    # p = (ln(n) + (eta - 1) * ln(ln(n))) / n
    log_n = math.log(n)
    log_log_n = math.log(log_n)

    p_numerator = log_n + (eta - 1) * log_log_n
    p = p_numerator / n
    
    # The final equation requires each number to be printed
    str_eta_minus_1 = f"({eta} - 1)"
    
    print("The d-threshold for Hamiltonicity is given by the formula:")
    print(f"p = (ln(n) + (η - 1) * ln(ln(n))) / n")
    print("\nFor the given values:")
    print(f"n = {n}")
    print(f"η = {eta}")
    
    print("\nThe equation is:")
    # Printing with the numbers plugged in
    print(f"p = (ln({n}) + {str_eta_minus_1} * ln(ln({n}))) / {n}")
    print(f"p = ({log_n:.4f} + {eta-1:.4f} * {log_log_n:.4f}) / {n}")
    print(f"p = ({log_n:.4f} + {(eta-1)*log_log_n:.4f}) / {n}")
    print(f"p = {p_numerator:.4f} / {n}")
    print(f"p = {p:.6f}")


# Example values for demonstration
n_example = 10000
eta_example = 20.0

calculate_hamiltonicity_threshold(n_example, eta_example)
<<<p = (math.log(10000) + (20.0 - 1) * math.log(math.log(10000))) / 10000>>>