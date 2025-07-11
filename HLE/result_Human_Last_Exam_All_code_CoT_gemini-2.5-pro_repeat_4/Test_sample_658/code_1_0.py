import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates the asymptotic value of the minimum degree d_B,delta
    for approximating e^-x on [0, B] with error delta.
    
    The asymptotic formula used is A(B, L) = B + L / log(L + 2),
    where L = log(1/delta).
    """
    if not (delta > 0 and delta < 1):
        raise ValueError("delta must be in the interval (0, 1)")
    if not (B >= 1):
        raise ValueError("B must be >= 1")

    # Calculate L
    L = math.log(1 / delta)

    # Calculate the asymptotic degree
    # Add a small epsilon to L in log to avoid math domain error if L is close to 0
    log_L_plus_2 = math.log(L + 2)
    asymptotic_d = B + L / log_L_plus_2

    # Print the equation with the specific numbers
    # We show the components of the calculation as requested.
    print(f"Given B = {B} and delta = {delta}")
    print(f"We calculate L = log(1/delta) = log(1/{delta}) = {L:.4f}")
    
    # Building the final equation string
    # B + L / log(L+2)
    equation_str = f"{B:.4f} + {L:.4f} / log({L:.4f} + 2)"
    result_str = f"= {B:.4f} + {L:.4f} / {log_L_plus_2:.4f} = {asymptotic_d:.4f}"
    
    print("\nThe asymptotic value is calculated using the formula: B + L/log(L+2)")
    print("Substituting the values:")
    print(equation_str)
    print(result_str)


# Example values for B and delta
B_val = 50
delta_val = 1e-10

calculate_asymptotic_degree(B_val, delta_val)