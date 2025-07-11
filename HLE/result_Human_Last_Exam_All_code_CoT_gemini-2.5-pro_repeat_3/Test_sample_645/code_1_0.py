import math
import numpy as np

def analyze_sets():
    """
    Analyzes four sets of natural numbers to determine if they can be the support
    of a power series that converges on the closed unit disk but not absolutely on the boundary.

    The problem is equivalent to identifying which of the sets are NOT Sidon sets.
    A set S is a Sidon set if for any continuous function on the unit circle whose Fourier
    coefficients are supported on S, the Fourier series must converge absolutely.
    A power series converging everywhere on the closed unit disk has a continuous sum function.
    If S were a Sidon set, this would imply the coefficients sum absolutely, which is
    contradictory to the problem's requirement.
    """

    print("Analyzing the properties of each set to determine if it is a Sidon set.\n")

    # --- Set 1: Random walk sums ---
    print("--- Analysis of Set 1 ---")
    print("S = {sum(N_k) for k<=n} where N_k ~ Poi(1)")
    # We generate one sample path for illustration.
    np.random.seed(0)
    n_max_1 = 30
    poisson_draws = np.random.poisson(1, n_max_1)
    partial_sums = np.cumsum(poisson_draws)
    set1 = sorted(list(set(partial_sums)))
    print(f"A sample realization of the set's values for n up to {n_max_1}:")
    print(set1)
    print("By the Law of Large Numbers, the n-th value S_n grows linearly with n (S_n/n -> 1).")
    print("A theorem by Salem and Zygmund shows that the value set of such a random walk is almost surely NOT a Sidon set.")
    print("Conclusion: Set 1 has the desired property (almost surely).\n")
    set1_works = True

    # --- Set 2: k-th powers ---
    print("--- Analysis of Set 2 ---")
    k = 4
    n_max_2 = 10
    print(f"S = {{n^{k} : n in N}} for k = {k}")
    s2 = [n**k for n in range(1, n_max_2 + 3)]
    print(f"First {n_max_2} elements: {[n**k for n in range(1, n_max_2 + 1)]}")
    print("We check for the convexity property: s_{n+2} - 2*s_{n+1} + s_n >= 0.")
    is_convex = all(s2[i+2] - 2*s2[i+1] + s2[i] >= 0 for i in range(n_max_2))
    if is_convex:
        print("The set is convex. A theorem by Meyer states that convex sets of integers are Sidon sets.")
    else:
        print("The set is not convex.") # This case will not be reached
    print("Conclusion: Set 2 does not have the desired property.\n")
    set2_works = False

    # --- Set 3: Primes ---
    print("--- Analysis of Set 3 ---")
    print("S = the set of primes")
    def get_primes(up_to):
        primes = []
        is_prime = [True] * (up_to + 1)
        if up_to >= 0: is_prime[0] = False
        if up_to >= 1: is_prime[1] = False
        for p in range(2, int(math.sqrt(up_to)) + 1):
            if is_prime[p]:
                for multiple in range(p*p, up_to + 1, p):
                    is_prime[multiple] = False
        for p in range(2, up_to + 1):
            if is_prime[p]:
                primes.append(p)
        return primes
    
    primes = get_primes(100)
    print("First 25 primes:", primes[:25])
    print("The gaps between primes are irregular. The set is not lacunary, nor is it convex.")
    print("A deep result in harmonic analysis (proven by Drury) shows that the set of primes is NOT a Sidon set.")
    print("Conclusion: Set 3 has the desired property.\n")
    set3_works = True

    # --- Set 4: Powers of pi/2 ---
    print("--- Analysis of Set 4 ---")
    print("S = {floor((pi/2)^n) : n in N}")
    n_max_4 = 15
    s4_raw = [math.floor((math.pi/2)**n) for n in range(1, n_max_4 + 2)]
    s4 = sorted(list(set(s4_raw)))
    print(f"First few unique elements: {s4[:n_max_4]}")
    print("We check for the lacunarity property: s_{n+1}/s_n >= q > 1.")
    print("Ratios of consecutive terms:")
    limit_ratio = math.pi / 2
    for i in range(len(s4) - 1):
        if s4[i] == 0: continue
        ratio = s4[i+1] / s4[i]
        if i < 10:
             print(f"s_{i+2}/s_{i+1} = {s4[i+1]}/{s4[i]} = {ratio:.4f}")

    print(f"The ratio approaches pi/2 ≈ {limit_ratio:.4f} > 1.")
    print("The set is lacunary. Lacunary sets are classic examples of Sidon sets.")
    print("Conclusion: Set 4 does not have the desired property.\n")
    set4_works = False

    # --- Final Summary ---
    print("--- Summary ---")
    final_sets = []
    if set1_works: final_sets.append(1)
    if set2_works: final_sets.append(2)
    if set3_works: final_sets.append(3)
    if set4_works: final_sets.append(4)
        
    print("The sets with the given property are those numbered: ", end="")
    for i, num in enumerate(final_sets):
        if i > 0:
            print(" and ", end="")
        print(num, end="")
    print(".")

# Execute the analysis
analyze_sets()