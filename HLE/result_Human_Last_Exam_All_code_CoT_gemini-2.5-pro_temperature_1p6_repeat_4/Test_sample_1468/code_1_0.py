import math

def calculate_lower_bound(N, q, d_prime):
    """
    Calculates the lower bound for the hidden layer dimension 'm' based on the derived formula.

    Args:
        N (int): The number of rows in the input matrix X.
        q (int): The sparsity parameter.
        d_prime (int): The dimension of the z_i vectors.

    Returns:
        int: The calculated lower bound for m.
    """
    # Problem constraints
    if not (q <= N / 2):
        print(f"Warning: Constraint q <= N/2 is not met ({q} > {N}/2). The derivation assumes this holds.")
    if not (q < d_prime):
        print(f"Warning: Constraint q < d' is not met ({q} >= {d_prime}). The derivation assumes this holds.")

    # The lower bound is derived from a dimensionality reduction argument.
    # The hidden layer must be large enough to encode k = floor(N/q) independent
    # projections, each of dimension d'.
    k = math.floor(N / q)
    lower_bound = k * d_prime
    
    # Print the equation and the result
    print("The derived lower bound for the hidden dimension m is given by the equation:")
    print(f"m >= floor(N / q) * d'")
    print("\nFor the given parameters:")
    print(f"N = {N}")
    print(f"q = {q}")
    print(f"d' = {d_prime}")
    print("\nThe calculation is:")
    print(f"m >= floor({N} / {q}) * {d_prime}")
    print(f"m >= {k} * {d_prime}")
    print(f"m >= {lower_bound}")
    print(f"\nThus, the asymptotic lower bound for m is Omega((N/q) * d').")

# Example values for the parameters
N_val = 1000
q_val = 10
d_prime_val = 50

# Execute the calculation and print the result
calculate_lower_bound(N_val, q_val, d_prime_val)
