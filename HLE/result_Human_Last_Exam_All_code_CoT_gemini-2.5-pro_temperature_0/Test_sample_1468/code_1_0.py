import math

def calculate_lower_bound(N, q, d_prime):
    """
    Calculates the lower bound for the hidden layer dimension m.

    Args:
        N (int): Number of input rows.
        q (int): Sparsity parameter.
        d_prime (int): Dimension of the z_i vectors.
    """
    # Based on the derivation, the lower bound for m is (N/q) * d'.
    # We assume N is a multiple of q for simplicity in the derivation,
    # but for a general lower bound, we use floor(N/q).
    M = math.floor(N / q)
    lower_bound = M * d_prime

    print(f"Given parameters:")
    print(f"N (number of rows) = {N}")
    print(f"q (sparsity) = {q}")
    print(f"d' (problem dimension) = {d_prime}")
    print("-" * 20)
    print(f"The lower bound for the hidden dimension 'm' is calculated as floor(N / q) * d'.")
    print(f"m >= floor({N} / {q}) * {d_prime}")
    print(f"m >= {M} * {d_prime}")
    print(f"m >= {lower_bound}")
    print("-" * 20)
    print(f"The asymptotic lower bound is Omega(N*d'/q).")
    return lower_bound

# Example values for the parameters, satisfying the problem constraints
# q <= N/2 and q < d'
N_val = 1000
q_val = 10
d_prime_val = 20

# Calculate and print the result
final_bound = calculate_lower_bound(N_val, q_val, d_prime_val)

# The question asks for the asymptotic lower bound.
# The derived bound is m >= (N/q) * d'.
# So the asymptotic bound is Omega(Nd'/q).
# The format <<<answer>>> requires a specific string.
# I will provide the formula as the answer.