import math

def calculate_minimum_curvature_cost(d, n):
    """
    Calculates and explains the minimum curvature cost for an NGD update.

    Args:
        d (int): The dimension of the square weight matrix (d x d).
        n (int): The number of training samples, where n < d.
    """
    if not n < d:
        print("Error: The condition n < d must be met.")
        return

    # The total number of parameters in the network is d*d.
    num_params = d**2
    
    # The matrix to be inverted in the NGD update is F + alpha*I,
    # where F is the Fisher Information Matrix. This matrix is (d^2 x d^2).
    # A naive inversion would cost O((d^2)^3) = O(d^6).
    
    # However, with n samples, the rank of the Fisher matrix F is at most n.
    # Since n < d, the matrix F is low-rank.
    
    # Using the Woodbury matrix identity, the inversion of the (d^2 x d^2) matrix
    # can be reduced to an inversion of an (n x n) matrix.
    
    # The cost of inverting an (n x n) matrix is O(n^3).
    # This is the minimum achievable cost for the inversion operation.
    
    # The final equation for the cost is n^3.
    base = n
    exponent = 3
    
    print("Problem Setup:")
    print(f"Layer dimension 'd': {d}")
    print(f"Number of samples 'n': {n} (note: n < d is satisfied)")
    print(f"Total number of parameters: {d}^2 = {num_params}")
    print("-" * 40)
    print("Analysis:")
    print(f"The NGD update requires inverting a {num_params}x{num_params} matrix.")
    print("By exploiting the low-rank structure (rank <= n), this can be simplified.")
    print("The minimum cost is achieved by inverting an n x n matrix instead.")
    print("-" * 40)
    print("Minimum Curvature Cost Calculation:")
    print("The complexity is O(n^3). The final equation for the cost is:")
    print(f"Cost is proportional to: {base}^{exponent}")
    
    # Calculate the proportional value for demonstration
    proportional_cost = base ** exponent
    print(f"Which evaluates to a value proportional to: {proportional_cost}")


# --- User-defined values ---
# d is the size of the layer (d x d)
d_value = 100
# n is the number of samples, with n < d
n_value = 30
# --- End of user-defined values ---

calculate_minimum_curvature_cost(d_value, n_value)