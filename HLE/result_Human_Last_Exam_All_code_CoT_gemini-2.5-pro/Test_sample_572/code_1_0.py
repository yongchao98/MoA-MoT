import math
import numpy as np

def main():
    """
    Main function to calculate r and demonstrate the FNP algorithm.
    """
    print("This script determines the largest rank 'r' for which an FNP algorithm can construct a rigid matrix.")
    print("------------------------------------------------------------------------------------------\n")

    # Part 1: Calculation of the largest r

    # We are given N and a small constant delta.
    # We will use example values for a concrete demonstration.
    N = 1000
    delta = 0.01
    # We'll assume the matrix is over the binary field GF(2), so the field size q=2.

    print("Step 1: The FNP Algorithm and Existence")
    print("An FNP algorithm can construct a rigid matrix by guessing a candidate and using an NP oracle to verify its rigidity.")
    print("This is possible if such matrices are known to exist.")
    print("A counting argument proves their existence for a sufficiently large r.")

    print("\nStep 2: Deriving the Bound for r")
    print("The counting argument shows that a matrix is rigid if r < N * (1 - sqrt(H_q(delta))).")
    print("Here, H_q(delta) is the q-ary entropy function.")
    print("For the binary field GF(2), q=2, the entropy H_2(delta) = -delta*log2(delta) - (1-delta)*log2(1-delta).")

    # Calculate H_2(delta)
    if delta <= 0 or delta >= 1:
        h_delta = 0
    else:
        h_delta = -delta * math.log2(delta) - (1 - delta) * math.log2(1 - delta)

    print(f"\nFor N = {N} and delta = {delta}:")
    print(f"H_2({delta}) = {h_delta:.4f}")

    # Calculate sqrt(H_2(delta))
    sqrt_h_delta = math.sqrt(h_delta)
    print(f"sqrt(H_2({delta})) = {sqrt_h_delta:.4f}")

    # Calculate the bound for r
    r_bound_factor = 1 - sqrt_h_delta
    r_bound = N * r_bound_factor
    r = int(math.floor(r_bound))

    print("\nStep 3: Calculating the Largest Rank r")
    print("The final equation for the largest rank r for which existence is known is:")
    # As requested, we output each number in the final equation
    print(f"r = floor(N * (1 - sqrt(H_2(delta))))")
    print(f"r = floor({N} * (1 - sqrt({h_delta:.4f})))")
    print(f"r = floor({N} * (1 - {sqrt_h_delta:.4f}))")
    print(f"r = floor({N} * {r_bound_factor:.4f})")
    print(f"r = floor({r_bound:.4f})")
    print(f"r = {r}")
    print("\nThis is the largest rank r for which this FNP construction method is known to work.")

    # Part 2: The FNP Algorithm (conceptual code)

    print("\n--- Illustrative FNP Algorithm Code ---")
    construct_rigid_matrix_fnp(N, delta, r)

# This is a hypothetical NP oracle. Its real-world implementation would be
# equivalent to solving an NP-complete problem.
def np_oracle_is_rank_lowerable(M, k, r):
    """
    Hypothetical NP oracle for the MATRIX_LOWER_RANK problem.
    Decides if matrix M's rank can be lowered to r by changing at most k entries.
    """
    print(f"\nOracle Call: Can rank of the given {M.shape[0]}x{M.shape[1]} matrix be lowered to {r} with <= {k} changes?")
    # The counting argument implies that for a randomly chosen matrix M,
    # the answer is 'False' with high probability for the r we are targeting.
    # We will simulate this probabilistic outcome.
    print("Oracle Response: False (the matrix is rigid, as expected for a random guess)")
    return False

def construct_rigid_matrix_fnp(N, delta, r):
    """
    Simulates the FNP algorithm to construct an (N x N) (delta, r)-rigid matrix.
    """
    print("\nExecuting FNP Algorithm:")
    # The non-deterministic "guess" is simulated by random generation.
    # A random matrix is a very good guess for a rigid matrix.
    print("1. Guessing a candidate matrix M (using randomness)...")
    # For display purposes, we'll create a smaller matrix if N is large.
    display_N = min(N, 8)
    M = np.random.randint(0, 2, size=(display_N, display_N))
    print("Candidate matrix M (a small version for display):")
    print(M)

    # The verification part using the NP oracle.
    print("\n2. Verifying the guess using an NP oracle...")
    k = int(delta * N**2)
    is_lowerable = np_oracle_is_rank_lowerable(M, k, r)

    if not is_lowerable:
        print("\n3. Verification successful! The algorithm halts and outputs the rigid matrix.")
        return M
    else:
        print("\n3. Verification failed. The guess was not a rigid matrix.")
        return None

if __name__ == '__main__':
    main()
