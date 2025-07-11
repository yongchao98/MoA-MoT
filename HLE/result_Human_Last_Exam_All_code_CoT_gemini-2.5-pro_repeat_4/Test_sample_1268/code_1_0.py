import math

def get_prime_factors(n):
    """
    Returns a set of prime factors of a natural number n.
    """
    factors = set()
    # Check for divisibility by 2
    while n % 2 == 0:
        factors.add(2)
        n = n // 2
    # Check for odd factors
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.add(i)
            n = n // i
    # If n is a prime number greater than 2
    if n > 2:
        factors.add(n)
    return factors

def is_squarefree(n):
    """
    Checks if a number is squarefree.
    """
    if n == 1:
        return True
    factors = get_prime_factors(n)
    for p in factors:
        if n % (p**2) == 0:
            return False
    return True

def calculate_norm_bound(N):
    """
    Calculates the upper bound for the maximum norm in relation to the covolume
    for a given squarefree natural number N.
    """
    print(f"Analyzing for the squarefree natural number N = {N}\n")

    if not isinstance(N, int) or N < 1:
        print("Error: N must be a natural number (integer >= 1).")
        return
        
    if not is_squarefree(N):
        print(f"Error: The number N = {N} is not squarefree.")
        return

    # 1. Get prime factors of N
    prime_factors = get_prime_factors(N)

    # 2. Calculate the index term psi(N) = N * product(1 + 1/p)
    psi_N = float(N)
    for p in prime_factors:
        psi_N *= (1.0 + 1.0/p)

    # 3. Calculate the covolume V
    V = (math.pi / 3.0) * psi_N

    # 4. Define the constant C = 9 / pi^2
    C = 9.0 / (math.pi**2)

    # 5. Calculate the final upper bound for the norm
    upper_bound = C * (V**2)

    # 6. Print the results step-by-step as requested
    print("Step 1: The covolume (V) is calculated based on N.")
    print(f"V = (π/3) * {N} * Π(1 + 1/p) for p|{N}")
    print(f"V = {V:.4f}\n")

    print("Step 2: The upper bound for the maximum norm (k_k,∞) is derived from V.")
    print("Derived relationship: k_k,∞ <= (9/π²) * V²\n")

    print("Final Calculation:")
    # The final equation with all numbers
    print(f"k_k,∞ <= {C:.4f} * ({V:.4f})²")
    print(f"k_k,∞ <= {upper_bound:.4f}")
    
    return upper_bound

# --- User Execution ---
# You can change the value of N to any squarefree natural number.
# Examples: 5, 6 (2*3), 10 (2*5), 30 (2*3*5), etc.
N_input = 30
final_bound = calculate_norm_bound(N_input)

# The final numerical answer is wrapped as requested
# The following line is for the final answer extraction.
if final_bound is not None:
    print(f"\n<<<k_k,∞ <= {final_bound:.4f}>>>")
