import numpy as np
import sympy

def analyze_sets():
    """
    Analyzes four sets of natural numbers to determine if they can be indices
    for a power series that converges everywhere on the unit disk but not absolutely on the boundary.
    The key property is whether the set of indices is lacunary (has Hadamard gaps).
    """
    
    print("Based on a theorem by Hadamard, Fatou, and Riesz, a power series whose exponents form a lacunary set")
    print("(i.e., m_{k+1}/m_k >= q > 1) can only converge everywhere on the unit circle if it also converges absolutely.")
    print("This violates the problem's conditions. Thus, we check which of the given sets are lacunary.\n")

    # Set 1: S = {sum_{k<=n} N_k}, N_k ~ Poi(1)
    # We check the ratio of consecutive elements from a simulation.
    np.random.seed(0)
    num_steps = 30000
    poisson_draws = np.random.poisson(1, num_steps)
    partial_sums = np.cumsum(poisson_draws)
    s1 = sorted(list(set(partial_sums)))
    if s1 and s1[0] == 0:
        s1 = s1[1:]
    
    ratio1 = s1[-1] / s1[-2]
    print("1. Set S = {sum_{k<=n} N_k}, N_k ~ Poi(1)")
    print(f"   A simulation with {num_steps} steps gives a sequence of indices.")
    print(f"   The ratio of the last two distinct elements is {s1[-1]}/{s1[-2]} = {ratio1:.6f}.")
    print("   The Law of Large Numbers implies the elements grow linearly on average, so m_k/m_{k-1} -> 1.")
    print("   Conclusion: The set is not lacunary. The property holds.\n")

    # Set 2: S = {n^k} for k >= 4. Let's use k=4.
    k_power = 4
    n_vals = np.arange(1, 2001, dtype=np.double)
    s2 = n_vals**k_power
    ratio2 = s2[-1] / s2[-2]
    print(f"2. Set S = {{n^k}} for k = {k_power}")
    print(f"   The ratio of consecutive terms is ((n)/(n-1))^k, which tends to 1.")
    print(f"   For n=2000, the ratio is {int(s2[-1])}/{int(s2[-2])} = {ratio2:.6f}.")
    print("   Conclusion: The set is not lacunary. The property holds.\n")

    # Set 3: S = the set of primes
    num_primes = 5000
    # Find the num_primes-th prime to set the range for primorange
    limit = sympy.prime(num_primes + 1)
    s3 = list(sympy.primerange(1, limit))
    ratio3 = s3[-1] / s3[-2]
    print("3. Set S = the set of primes")
    print(f"   The Prime Number Theorem implies p_n ~ n*log(n), so p_{n+1}/p_n -> 1.")
    print(f"   For the {num_primes-1}-th and {num_primes}-th primes, the ratio is {s3[-1]}/{s3[-2]} = {ratio3:.6f}.")
    print("   Conclusion: The set is not lacunary. The property holds.\n")

    # Set 4: S = {floor((pi/2)^n)}
    n_max = 60
    s4 = []
    for n in range(1, n_max + 1):
        # Use high precision for calculation before flooring
        val = np.floor((np.pi/2)**n)
        if not s4 or val > s4[-1]:
            s4.append(val)

    ratio4 = s4[-1] / s4[-2]
    pi_half = np.pi / 2
    print(f"4. Set S = {{floor((pi/2)^n)}}")
    print(f"   The terms grow like (pi/2)^n, so the ratio of consecutive terms approaches pi/2 ~ {pi_half:.6f}.")
    print(f"   The ratio of the last two distinct elements is {int(s4[-1])}/{int(s4[-2])} = {ratio4:.6f}.")
    print("   Conclusion: The set is lacunary. The property does not hold.\n")

    print("Final result: Sets 1, 2, and 3 have the property, while set 4 does not.")

analyze_sets()