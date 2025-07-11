import math

def demonstrate_query_lower_bound():
    """
    Explains and demonstrates the super-polynomial lower bound for learning
    two-hidden-layer ReLU networks in the SQ model.
    """
    # The theoretical lower bound for the number of queries (Q) is d^Ω(log d).
    # This means there exists some positive constant 'c' such that for large d,
    # the number of queries is at least d^(c * log(d)).
    # For this demonstration, we will use c=1 and the base-2 logarithm.
    
    c = 1

    print("The minimum number of queries (Q) required to learn the specified neural network is super-polynomial in the dimension d.")
    print("This complexity is expressed by the lower bound: Q(d) = d^Ω(log d).")
    print("\nTo illustrate, we can write out the equation for a specific instance of this bound.")
    print("Let's use the constant c=1 and the base-2 logarithm for our example equation.")
    
    # As requested, printing each number in the final equation. The only number is the constant 'c'.
    print("\n--- Demonstrative Equation ---")
    print(f"Q(d) >= d^({c} * log₂(d))")
    print("----------------------------")

    print("\nThis value grows extremely fast. Here are some examples:")

    # A helper function to calculate the lower bound for a given d
    def calculate_bound(d_val, c_val):
        if d_val <= 1:
            return 1
        # Use log base 2 for the calculation
        log2_d = math.log2(d_val)
        # The exponent is c * log2(d)
        exponent = c_val * log2_d
        # The result is d^(exponent), which can be calculated as 2^(log2(d) * exponent)
        # to maintain precision and avoid potential overflow with large bases.
        return 2**(log2_d * exponent)

    # Demonstrate the growth for a few values of d
    for d in [8, 16, 32, 64, 128]:
        # Calculate the lower bound value
        query_count = calculate_bound(d, c)
        print(f"For d = {d:3}, the minimum query count is on the order of {query_count:,.2e}")

# Execute the demonstration
demonstrate_query_lower_bound()