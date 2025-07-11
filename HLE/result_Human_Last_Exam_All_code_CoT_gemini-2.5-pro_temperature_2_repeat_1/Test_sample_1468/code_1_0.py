import math

def calculate_lower_bound(N, d_prime, q):
    """
    Calculates the lower bound for m based on the derived formula.

    Args:
        N (int): Total number of data points.
        d_prime (int): Dimension of the feature vectors z_i.
        q (int): Sparsity of the average.
    """
    # --- Sanity checks based on problem description ---
    valid = True
    if not (q <= N / 2):
        print(f"Warning: The condition q <= N/2 is not met. (q={q}, N/2={N/2})")
        valid = False
    if not (q < d_prime):
        print(f"Warning: The condition q < d' is not met. (q={q}, d'={d_prime})")
        valid = False
    if q < 1:
        print(f"Warning: q must be >= 1. (q={q})")
        valid = False
    if not valid:
        print("The provided parameters may not be valid for this problem.\n")

    # --- Calculation ---
    # The lower bound for m is given by the formula: m >= d' * floor(N / (q + 1))
    
    # Check if q+1 > N, which would make the floor result 0
    if N < (q + 1):
        k_val = 0
    else:
        k_val = math.floor(N / (q + 1))
    
    lower_bound = d_prime * k_val

    # --- Output ---
    print("Derivation of the lower bound for m:")
    print("="*40)
    print("Problem parameters:")
    print(f"  N (number of points) = {N}")
    print(f"  q (sparsity) = {q}")
    print(f"  d' (feature dimension) = {d_prime}")
    print("")
    print("The lower bound is calculated using the formula:")
    print("  m >= d' * floor(N / (q + 1))")
    print("")
    print("Plugging in the numbers:")
    print(f"  m >= {d_prime} * floor({N} / ({q} + 1))")
    print(f"  m >= {d_prime} * floor({N} / {q+1})")
    print(f"  m >= {d_prime} * {k_val}")
    print(f"  m >= {lower_bound}")
    print("="*40)
    print(f"The calculated lower bound for m is: {lower_bound}")
    print(f"Asymptotically, this corresponds to m = Omega(N * d' / q).")


if __name__ == '__main__':
    # --- User-defined parameters for an example calculation ---
    N = 1000
    q = 10
    d_prime = 50
    
    calculate_lower_bound(N, d_prime, q)