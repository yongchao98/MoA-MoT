import math

def h2_inv_half():
    """
    Solves H2(x) = 0.5 for x in [0, 0.5]. This value is a constant.
    H2(x) = -x*log2(x) - (1-x)*log2(1-x)
    A simple numerical solver (binary search).
    """
    y = 0.5
    low = 0.0
    high = 0.5
    for _ in range(100): # 100 iterations for high precision
        mid = (low + high) / 2
        if mid == 0:
            h = 0
        else:
            h = -mid * math.log2(mid) - (1 - mid) * math.log2(1 - mid)
        
        if h < y:
            low = mid
        else:
            high = mid
    return low

def calculate_rigidity_rank(N, delta):
    """
    Calculates the maximum rank 'r' for which we can construct a (delta, r)-rigid matrix
    based on the method of building a matrix from a good code.
    
    We use a binary code with a rate of 1/2.
    """
    # 1. We use FNP to construct a code with rate 1/2.
    code_rate = 0.5
    K = math.floor(code_rate * N)
    
    # 2. From Gilbert-Varshamov bound, for a code with rate 0.5, we can achieve a
    # relative distance `beta` where H2(beta) <= 1 - 0.5 = 0.5.
    # We find beta for H2(beta) = 0.5.
    beta = h2_inv_half()
    D = math.floor(beta * N)
    
    # 3. An attacker wants to reduce the rank from K to r. This requires at least (K-r)D changes.
    # For the matrix to be (delta, r)-rigid, this must be more than delta*N^2.
    # (K - r) * D > delta * N^2
    # K - r > (delta * N^2) / D
    # r < K - (delta * N^2) / D
    # r < K - (delta * N^2) / (beta * N) = K - (delta * N) / beta
    
    r_bound = K - (delta * N) / beta
    if r_bound < 0:
        r_bound = 0
        
    print("This script calculates the largest rank 'r' for which the described FNP construction yields a (delta, r)-rigid matrix.")
    print("-" * 50)
    print("Input parameters for the rigidity problem:")
    print(f"N (matrix size) = {N}")
    print(f"delta (fraction of allowed changes) = {delta}")
    print("-" * 50)
    
    print("Parameters of the FNP-constructed matrix:")
    print(f"Underlying code has rate (K/N) = {code_rate}")
    print(f"The resulting matrix rank is K = {K}")
    print(f"The relative distance of the code is beta, where H_2(beta) = {1.0 - code_rate}. Thus beta ~= {beta:.4f}")
    print(f"The absolute minimum distance of the code is D ~= {D}")
    print("-" * 50)

    print("The final equation for the rank 'r' is derived from the condition:")
    print("Number of changes to reduce rank > delta * N^2")
    print("(K - r) * D > delta * N^2")
    print("\nWhich simplifies to the bound: r < K - (delta * N) / beta")
    print("Plugging in the numbers:")
    num_N = N
    num_K = K
    num_delta = delta
    num_beta = beta
    
    print(f"r < {num_K} - ({num_delta} * {num_N}) / {num_beta:.4f}")
    print(f"r < {r_bound:.2f}")
    
    final_r = math.floor(r_bound)
    print(f"\nThus, for the given N and delta, we can construct a ({delta}, {final_r})-rigid matrix.")
    
# --- Main Execution ---
# You can change these values to see how the result changes.
N = 10000
delta = 0.005

calculate_rigidity_rank(N, delta)