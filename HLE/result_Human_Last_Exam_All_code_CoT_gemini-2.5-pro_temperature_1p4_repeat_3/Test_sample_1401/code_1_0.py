import math

def illustrate_sq_lower_bound(d, c=1.0):
    """
    Illustrates the quasi-polynomial lower bound for learning two-layer ReLU networks
    using Statistical Query (SQ) algorithms.

    Args:
        d (int): The input dimension of the problem.
        c (float): A constant factor representing the hidden constants in the
                 Omega notation of the complexity bound.
    """
    # Introduction to the problem and the theoretical result
    print("Problem: Minimum number of Statistical Queries (SQ) to learn a poly(d)-sized,")
    print("two-hidden-layer ReLU network on N(0, I_d) inputs to a squared loss of 1/poly(d).")
    print("-" * 70)
    print("From theoretical computer science, the lower bound on the number of queries for this problem is known to be quasi-polynomial.")
    print("A quasi-polynomial function grows faster than any polynomial but slower than an exponential function.")
    print("The specific form of this lower bound is d^Ω(log d), which we can write as d^(c * log(d)) for some constant c > 0.")
    print("-" * 70)

    # Check for valid input dimension
    if d <= 1:
        print("The dimension 'd' must be greater than 1 for a meaningful bound.")
        return

    # Calculation
    # We use log base 2 for convention, as is common in complexity theory.
    log2_d = math.log2(d)
    exponent_value = c * log2_d
    
    # Calculate the lower bound value.
    # The formula is d^exponent_value. To compute this, we can use powers or
    # the identity x^y = exp(y * ln(x)).
    log_e_d = math.log(d)
    lower_bound_value = math.exp(exponent_value * log_e_d)

    # Output the breakdown of the final equation
    print(f"For dimension d = {d} and a sample constant c = {c}:\n")
    print("The final equation for the lower bound has the form: base ^ (c * log2(base))")
    print(f"    base          = {d}")
    print(f"    c             = {c}")
    print(f"    log2(base)    = log2({d}) ≈ {log2_d:.4f}")
    print(f"    exponent      = {c} * {log2_d:.4f} ≈ {exponent_value:.4f}")
    print("-" * 70)
    print(f"Final Calculation: {d}^({exponent_value:.4f})")
    print(f"Minimum Number of Queries >= {lower_bound_value:e}")


# --- Main execution ---
# You can change this value to see how the lower bound grows with the dimension.
input_dimension = 10
illustrate_sq_lower_bound(input_dimension)