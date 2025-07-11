import math

def analyze_quadratic_field(N):
    """
    Calculates and demonstrates the relationship between the Minkowski bound and
    the covolume for a real quadratic field Q(sqrt(N)).
    """
    if N <= 0 or int(math.sqrt(N)) ** 2 == N:
        print(f"Error: N must be a non-squarefree positive integer. N={N} is not valid.")
        # Check for squarefree property
        for i in range(2, int(math.sqrt(N)) + 1):
            if N % (i**2) == 0:
                print(f"Error: N must be squarefree. N={N} is divisible by {i**2}.")
                return
        return
        
    print(f"Analyzing the real quadratic field K = Q(sqrt({N}))\n")

    # 1. Determine parameters
    n = 2  # degree
    r2 = 0 # number of complex embedding pairs

    # 2. Calculate the discriminant
    if N % 4 == 1:
        delta_K = N
    else: # N % 4 is 2 or 3
        delta_K = 4 * N
    
    print(f"The discriminant is Delta_K = {delta_K}")

    # 3. Calculate the covolume (V)
    # V = 2^(-r2) * sqrt(|Delta_K|) = 2^0 * sqrt(|Delta_K|)
    V = math.sqrt(delta_K)
    print(f"The covolume is V = sqrt({delta_K}) \u2248 {V:.4f}")

    # 4. Calculate the Minkowski Bound (M_K)
    # M_K = (4/pi)^r2 * n!/n^n * sqrt(|Delta_K|)
    # M_K = (4/pi)^0 * 2!/2^2 * sqrt(|Delta_K|) = (1/2) * sqrt(|Delta_K|)
    M_K = (1 / 2) * math.sqrt(delta_K)
    print(f"The Minkowski bound is M_K = (1/2) * sqrt({delta_K}) \u2248 {M_K:.4f}")

    # 5. Show the relationship k_k_inf <= (1/2) * V
    print("\n" + "="*50)
    print("The theoretical upper bound for the maximum norm (k_k_inf) is M_K.")
    print("The derived relationship is: k_k_inf <= (1/2) * V")
    
    # Show the final equation with calculated numbers
    coefficient = 1/2
    print("\nDemonstrating with the calculated values:")
    print(f"M_K <= {coefficient} * V")
    print(f"{M_K:.4f} <= {coefficient} * {V:.4f}")
    print(f"{M_K:.4f} <= {coefficient * V:.4f} (This holds true)")

    print("\n" + "="*50)
    print("Final Answer: The upper bound is expressed by the equation:")
    print(f"k_k_inf <= ({1} / {2}) * V")
    

# Example with a squarefree natural number, N=10 (10 = 2*5, and 10 % 4 == 2)
# You can change N to any other squarefree positive integer like 3, 5, 7, 13, 15 etc.
N = 10
analyze_quadratic_field(N)
