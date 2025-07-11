import math

def calculate_asymptotic_degree(B, delta):
    """
    Calculates and explains the asymptotic value A(B, delta) for d_{B,delta}.
    The formula is A(B, delta) = B + L / log(1 + L/B), where L = log(1/delta).
    This script demonstrates the calculation with the given B and delta values.
    """
    print(f"This script calculates the asymptotic value of d_{B,delta} for B = {B} and delta = {delta}.")
    
    if not (B >= 1 and 0 < delta < 1):
        print("Error: Input values are out of the specified range (B >= 1, 0 < delta < 1).")
        return

    # Calculate L = log(delta^{-1})
    L = math.log(1 / delta)
    
    # The final formula is A = B + L / log(1 + L/B)
    # The script will now print the components of this formula step-by-step.
    
    print("\nThe final equation is of the form: Result = B + L / log(1 + L/B)")
    print("Let's substitute the numbers into the equation:")
    
    # Calculate the term inside the logarithm
    log_term_inner = 1 + L / B
    
    # Calculate the denominator of the fraction
    log_term_outer = math.log(log_term_inner)
    
    # Calculate the fractional part of the expression
    # A check is included for the case where the denominator is zero
    if log_term_outer == 0:
        # This case occurs when L/B is extremely small.
        # We use the approximation log(1+x) ≈ x for small x.
        # The fractional part L / log(1+L/B) becomes L / (L/B) = B.
        fractional_part = B
    else:
        fractional_part = L / log_term_outer
    
    # Calculate the final result
    result = B + fractional_part
    
    print(f"Result = {B:.4f} + {L:.4f} / log(1 + {L:.4f}/{B:.4f})")
    print(f"Result = {B:.4f} + {L:.4f} / log({log_term_inner:.4f})")
    print(f"Result = {B:.4f} + {L:.4f} / {log_term_outer:.4f}")
    print(f"Result = {B:.4f} + {fractional_part:.4f}")
    print(f"Result = {result:.4f}")


# Using example values to demonstrate the calculation.
# Let's choose values that are in the crossover regime, for instance B=10 and L=10.
B_val = 10.0
# L = 10 implies delta = e^{-10}
delta_val = math.exp(-10)

calculate_asymptotic_degree(B=B_val, delta=delta_val)