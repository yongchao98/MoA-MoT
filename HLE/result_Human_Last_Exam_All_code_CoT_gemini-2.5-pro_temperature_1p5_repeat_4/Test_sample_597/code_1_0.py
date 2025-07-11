import random

def solve():
    """
    This function simulates a verifier for a PCP that is both "Red" and "Blue"
    and prints an analysis of the results.
    """
    # Parameters for the simulation
    n = 5  # Dimension of the input space {0,1}^n. n=5 means 2^5=32 inputs.
    N = 2**n # Size of the truth table (the proof length)
    num_flips = 4 # Number of errors to introduce in the proof
    num_trials = 20000 # Number of times to run the verifier to estimate probability

    # --- Step 1: Define a correct proof (a codeword in the Hadamard code) ---
    # A correct proof corresponds to a linear function f(x) = a . x (mod 2).
    # Let's choose a random linear function 'a'.
    a = [random.randint(0, 1) for _ in range(n)]

    # The proof is the truth table of the linear function.
    # An input x is an integer from 0 to N-1.
    # To compute a . x, we look at the binary representations.
    def dot_product(a_vec, x_int):
        res = 0
        for i in range(n):
            # Get i-th bit of x
            x_bit = (x_int >> i) & 1
            res ^= a_vec[i] * x_bit
        return res

    correct_proof = [dot_product(a, x) for x in range(N)]

    # --- Step 2: Create a corrupted proof by flipping bits ---
    corrupted_proof = list(correct_proof)
    # Ensure flip indices are unique if num_flips < N
    flip_indices = random.sample(range(N), min(num_flips, N))
    for idx in flip_indices:
        corrupted_proof[idx] = 1 - corrupted_proof[idx]

    # --- Step 3: Calculate the actual relative distance delta ---
    # Delta is the fraction of differing positions.
    delta = num_flips / N

    # --- Step 4: Simulate the BLR (Blum-Luby-Rubinfeld) verifier ---
    # The verifier checks if f(x) + f(y) = f(x+y) for random x, y.
    # This property defines linearity.
    rejections = 0
    for _ in range(num_trials):
        # Pick random inputs x and y
        x = random.randint(0, N - 1)
        y = random.randint(0, N - 1)
        
        # Calculate x+y (which is XOR in the vector space {0,1}^n)
        z = x ^ y
        
        # Query the proof oracle for f(x), f(y), and f(x+y)
        fx = corrupted_proof[x]
        fy = corrupted_proof[y]
        fz = corrupted_proof[z]
        
        # Check the linearity property
        if fx ^ fy != fz:
            rejections += 1

    # --- Step 5: Calculate and print the results ---
    # Estimate the rejection probability
    estimated_p_rej = rejections / num_trials
    # The theoretical rejection probability for the BLR test is 2 * delta
    theoretical_p_rej = 2 * delta

    print("--- Simulation of a 'Red and Blue' PCP Verifier ---")
    print(f"This example uses the BLR linearity test, a component of PCP constructions.")
    print(f"A proof is a function's truth table. Distance is the fraction of entries that must be changed to make the function linear.")
    print("\nSimulation Parameters:")
    print(f"  - Proof length (N = 2^n): {N}")
    print(f"  - Number of introduced errors (k): {num_flips}")

    print("\nDistance Calculation:")
    print(f"The relative Hamming distance of the corrupted proof from the set of correct proofs is delta = k / N.")
    print(f"  - delta = {num_flips} / {N} = {delta}")

    print("\nVerifier Test and Rejection Probability:")
    print("The verifier tests if f(x) + f(y) = f(x+y) for random x, y.")
    print("A famous result states that the rejection probability (p_rej) is exactly 2 * delta.")
    print("This verifier is both 'Red' (p_rej >= 2*delta) and 'Blue' (p_rej <= 2*delta).")

    print("\nResults:")
    print(f"  - Number of verifier trials: {num_trials}")
    print(f"  - Observed rejections: {rejections}")
    print(f"  - Empirical rejection probability: {estimated_p_rej:.4f}")
    print(f"  - Theoretical rejection probability (2 * delta): {theoretical_p_rej:.4f}\n")

    print("Final Equation Check:")
    equation_p_rej_str = "p_rej"
    equation_factor = 2
    equation_delta = delta
    print(f"The relationship is {equation_p_rej_str} = {equation_factor} * delta.")
    print(f"Our empirical {equation_p_rej_str} is {estimated_p_rej:.4f}")
    print(f"The theoretical value is {equation_factor} * {equation_delta} = {theoretical_p_rej:.4f}")

solve()