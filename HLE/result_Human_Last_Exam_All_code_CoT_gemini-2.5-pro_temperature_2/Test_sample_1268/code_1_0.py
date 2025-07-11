import math

def get_minkowski_relation(N):
    """
    Calculates and demonstrates the relationship between the Minkowski bound (k_k,inf)
    and the covolume (V) for a real quadratic field Q(sqrt(N)).

    Args:
        N (int): A squarefree natural number.
    """
    print(f"Analyzing the field Q(sqrt({N}))...\n")

    # Step 1: Calculate the discriminant (delta_k)
    if N % 4 == 1:
        delta_k = N
    else: # N % 4 is 2 or 3
        delta_k = 4 * N
    
    print(f"The discriminant (Delta_K) is: {delta_k}")

    # Step 2: Calculate the covolume (V)
    # For a real quadratic field, V = sqrt(Delta_K)
    covolume_V = math.sqrt(delta_k)
    print(f"The covolume (V) is sqrt({delta_k}) = {covolume_V:.4f}")

    # Step 3: Calculate the Minkowski bound (k_k,inf)
    # For a real quadratic field, k_k,inf = (1/2) * sqrt(Delta_K)
    k_k_inf = 0.5 * math.sqrt(delta_k)
    print(f"The Minkowski bound (k_k,inf) is (1/2) * sqrt({delta_k}) = {k_k_inf:.4f}\n")

    # Step 4: Display the final relationship
    # The relationship is k_k,inf = 0.5 * V
    coefficient = 0.5
    print("The relationship between the upper bound (k_k,inf) and the covolume (V) is k_k,inf = 0.5 * V.")
    print("For our example:")
    
    # Final print of the equation with the calculated numbers
    print(f"{k_k_inf:.4f} = {coefficient} * {covolume_V:.4f}")

# Example usage with N=5 (a squarefree natural number)
# 5 mod 4 = 1, so Delta_K = 5
# V = sqrt(5) approx 2.236
# k_k,inf = 0.5 * sqrt(5) approx 1.118
get_minkowski_relation(5)
