import math

def calculate_bound_for_max_norm(N):
    """
    Calculates the upper bound for the maximum norm (k_k,inf) in relation
    to the covolume (V) for the ring of integers of the real quadratic
    field Q(sqrt(N)).

    Args:
        N (int): A squarefree natural number.
    """
    # --- Step 1: Explain the setup ---
    print(f"Analyzing the case for N = {N}.")
    print(f"This corresponds to the real quadratic number field K = Q(sqrt({N})).")
    print("We consider the lattice formed by the ring of integers of this field.")
    print("The maximum norm, k_k,inf, is interpreted as the L-infinity norm of the shortest non-zero vector in this lattice.")
    print("The covolume is denoted by V.")
    print("-" * 30)

    # --- Step 2: Calculate the Discriminant (Delta) ---
    if N % 4 == 1:
        delta = N
        print(f"Since N % 4 == 1, the discriminant is Delta = N = {delta}.")
    else: # N % 4 == 2 or 3
        delta = 4 * N
        print(f"Since N % 4 != 1, the discriminant is Delta = 4 * N = {delta}.")

    # --- Step 3: Calculate the Covolume (V) ---
    V = math.sqrt(delta)
    print(f"The covolume of the lattice is V = sqrt(Delta) = sqrt({delta}) = {V:.4f}.")

    # --- Step 4: Calculate the Upper Bound ---
    # The relationship from Minkowski's theorem is k_k,inf <= sqrt(V)
    upper_bound = math.sqrt(V)
    print(f"The theoretical upper bound for k_k,inf is sqrt(V).")
    print("-" * 30)

    # --- Step 5: Output the final equation with all numbers ---
    print("Final Inequality:")
    # Using f-string to format the final equation
    print(f"k_k,inf <= sqrt(V)")
    print(f"k_k,inf <= sqrt({V:.4f})")
    print(f"k_k,inf <= {upper_bound:.4f}")
    
# Let's use N=7 as an example, which is a squarefree natural number.
calculate_bound_for_max_norm(7)