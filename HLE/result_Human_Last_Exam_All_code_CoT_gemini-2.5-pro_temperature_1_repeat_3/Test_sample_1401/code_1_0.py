import math

def calculate_exponential_lower_bound(d, c=0.1):
    """
    Calculates an illustrative value for the exponential lower bound on SQ queries.

    The theoretical lower bound for this learning problem is exp(Omega(d)),
    meaning it grows exponentially with the dimension d. We use a specific
    function, exp(c * d), to demonstrate this behavior with a concrete number.

    Args:
        d (int): The input dimension of the neural network.
        c (float): A constant representing the hidden factor in the Omega (Ω) notation.
                   A smaller c represents a more optimistic (but still exponential) bound.

    Returns:
        float: The calculated lower bound on the number of queries.
    """
    # The minimum number of queries is proven to be exponential in dimension d.
    min_queries_value = math.exp(c * d)

    print("The theoretical lower bound on the number of queries is of the form: exp(c * d)")
    print("For the chosen parameters, the equation is:")
    # The final print statement shows each number in the equation, as requested.
    print(f"exp({c} * {d}) = {min_queries_value}")

    return min_queries_value


if __name__ == '__main__':
    # --- Problem Setup ---
    # We consider a moderately high dimension to see the exponential growth.
    input_dimension = 100

    print(f"Problem: Find the minimum number of SQ queries to learn a two-hidden-layer ReLU network.")
    print(f"Input Dimension (d): {input_dimension}")
    print("Condition: Query tolerance is not negligible in d (e.g., > 1/poly(d)).\n")

    print("--- Calculation ---")
    # Calculate and print the illustrative lower bound.
    # The constant 'c' in exp(c*d) is chosen as 0.1 for this example.
    calculate_exponential_lower_bound(d=input_dimension, c=0.1)