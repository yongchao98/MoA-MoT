import numpy as np
import math

def analyze_sets():
    """
    Investigates the properties of the four sets to determine if they are
    finite unions of arithmetic progressions, which is the criterion for
    solving the problem.
    """
    limit_n = 50
    limit_x = 1000000

    print("Analyzing the properties of each set:")
    print("-" * 40)

    # --- Case 1: S = {sum(N_k) for k<=n}, N_k ~ Poi(1) ---
    # A finite union of APs has a finite set of differences between consecutive elements.
    # We check if the set of differences for a realization of S is finite.
    print("Case 1: S = {sum(N_k) for k<=n}, N_k ~ Poi(1)")
    np.random.seed(0)  # for reproducibility
    poisson_draws = np.random.poisson(1, size=limit_n)
    s1 = np.cumsum(poisson_draws)
    s1 = sorted(list(set(s1))) # remove duplicates, as Poi(1) can be 0
    diffs1 = np.diff(s1)
    print(f"A realization of the first few elements of S: {s1[:15]}")
    print(f"Differences between these elements: {diffs1[:14]}")
    print(f"The set of unique differences found is {sorted(list(set(diffs1)))}.")
    print("Since the differences are drawn from a Poisson distribution, the set of all differences is almost surely infinite. Thus, S is not a finite union of APs.\n")

    # --- Case 2: S = {n^k} for k=4 ---
    # A finite union of APs has a positive asymptotic density.
    # We check if the density of S approaches 0.
    print("Case 2: S = {n^k} for k=4")
    k = 4
    s2_count_at_x = lambda x: math.floor(x**(1/k))
    for x_power in [2, 4, 6]:
        x = 10**x_power
        count = s2_count_at_x(x)
        density = count / x
        print(f"For X = {x:10d}, |S intersect [1,X]| = {count:5d}. Density = {count}/{x} = {density:.6f}")
    print("The density, proportional to X^(1/k-1), approaches 0. Thus, S is not a finite union of APs.\n")

    # --- Case 3: S = the set of primes ---
    # We check if the density of primes approaches 0.
    print("Case 3: S = the set of primes")
    def count_primes(n):
        if n < 2: return 0
        primes = [True] * (n + 1)
        primes[0] = primes[1] = False
        for i in range(2, int(math.sqrt(n)) + 1):
            if primes[i]:
                for multiple in range(i*i, n + 1, i):
                    primes[multiple] = False
        return sum(primes)

    for x_power in [2, 4, 6]:
        x = 10**x_power
        count = count_primes(x)
        density = count / x
        print(f"For X = {x:10d}, pi(X) = {count:5d}. Density = {count}/{x} = {density:.6f}")
    print("The density, proportional to 1/ln(X), approaches 0. Thus, S is not a finite union of APs.\n")

    # --- Case 4: S = {floor((pi/2)^n)} ---
    # This is a lacunary set. We check its density.
    print("Case 4: S = {floor((pi/2)^n)}")
    q = math.pi / 2
    s4_count_at_x = lambda x: math.floor(math.log(x, q)) if x > 1 else 0
    for x_power in [2, 4, 6]:
        x = 10**x_power
        count = s4_count_at_x(x)
        density = count / x
        print(f"For X = {x:10d}, |S intersect [1,X]| = {count:5d}. Density = {count}/{x} = {density:.6f}")
    print("The density, proportional to log(X)/X, approaches 0. Thus, S is not a finite union of APs.\n")

    print("-" * 40)
    print("Conclusion: All four sets are not finite unions of arithmetic progressions.")
    print("By the Rudin-Salem theorem, all four sets have the desired property.")

# Execute the analysis
analyze_sets()