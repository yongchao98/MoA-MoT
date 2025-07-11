import math

def illustrate_sq_lower_bound(d, c=1.0):
    """
    Illustrates the super-polynomial SQ lower bound for learning certain
    two-hidden-layer ReLU networks.

    The theoretical lower bound on the number of queries is d^Ω(log d).
    This function calculates d^(c * log(d)) for a specific constant c
    to demonstrate this rapid growth.

    Args:
        d (int): The dimension of the input space.
        c (float): A constant factor for the exponent, representing the
                   hidden constant in the big-Omega Ω notation.
    """
    if d <= 1:
        print(f"Dimension d must be greater than 1. Received: {d}")
        return

    # In theoretical hardness reductions, the learning problem is often shown
    # to be as hard as learning a parity function on k variables, where k
    # is chosen to be logarithmic in the dimension d.
    # We choose k = c * log2(d).
    try:
        k = c * math.log2(d)
    except ValueError:
        print(f"Invalid value for d: {d}. Must be positive.")
        return

    # The SQ lower bound for learning the corresponding hard function is d^Ω(k).
    # We substitute k = c * log2(d) to get d^Ω(log d).
    # The calculation below shows this value for our choice of d.
    try:
        lower_bound_queries = math.pow(d, k)
    except OverflowError:
        # The number can grow extremely fast, quickly exceeding standard float limits.
        lower_bound_queries = float('inf')

    print(f"--- Illustration for dimension d = {d} ---")
    print("The minimum number of queries is given by the formula: d^Ω(log d)")
    print("This means the number of queries grows faster than any polynomial in d.")
    print("\nTo demonstrate, we calculate d^(c * log2(d)) for c=1.0:")
    
    # Per the instructions, we output each number in the final equation.
    print("\n[Final Equation Numbers]")
    print(f"Base (d): {d}")
    print(f"Exponent (k = c * log2(d)): {k:.4f}")
    print(f"Result (d^k): {lower_bound_queries:e}")
    
    print(f"\nEquation: {d:.1f} ^ {k:.4f} = {lower_bound_queries:e}\n")


# --- Main execution ---
# Set an example dimension d to see the complexity.
# Notice how quickly the number of queries grows.
example_dimension = 10
illustrate_sq_lower_bound(example_dimension)

# A slightly larger dimension
example_dimension_2 = 30
illustrate_sq_lower_bound(example_dimension_2)
