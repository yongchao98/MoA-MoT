import numpy as np
import math

def get_set_1(n_terms=20):
    """
    Generates terms for S = {sum_{k<=n} N_k} where N_k ~ Poi(1).
    This is a random realization of the set.
    """
    # N_k are i.i.d. Poi(1)
    N_k = np.random.poisson(1, size=n_terms)
    s = np.cumsum(N_k)
    # Ensure the set has unique, sorted elements
    s = sorted(list(set(s)))
    return s

def get_set_2(n_terms=20, k=4):
    """
    Generates terms for S = {n^k} for k>=4.
    """
    s = [n**k for n in range(1, n_terms + 1)]
    return s

def get_set_3(n_terms=20):
    """
    Generates the first n_terms primes.
    """
    primes = []
    num = 2
    # Simple prime checker
    def is_prime(m):
        if m < 2:
            return False
        for i in range(2, int(math.sqrt(m)) + 1):
            if m % i == 0:
                return False
        return True
    
    while len(primes) < n_terms:
        if is_prime(num):
            primes.append(num)
        num += 1
    return primes

def get_set_4(n_terms=20):
    """
    Generates terms for S = {floor((pi/2)^n)}.
    """
    q = math.pi / 2
    s = [math.floor(q**n) for n in range(1, n_terms + 1)]
    # Ensure the set has unique, sorted elements as floor can yield duplicates
    s = sorted(list(set(s)))
    return s

def compute_gaps(s):
    """
    Computes the gaps between consecutive elements of a list.
    """
    if len(s) < 2:
        return []
    return [s[i+1] - s[i] for i in range(len(s) - 1)]

def main():
    """
    Main function to analyze the four sets.
    """
    print("Analyzing the four sets to check for the 'unbounded gaps' property.")
    print("A sufficient condition for the existence of the described power series is that the set of indices S has unbounded gaps.\n")

    # Set 1
    print("--- Set 1: S = {sum_{k<=n} N_k}, N_k ~ Poi(1) ---")
    s1 = get_set_1(20)
    gaps1 = compute_gaps(s1)
    print(f"A random realization of the first {len(s1)} terms:\n{s1}")
    print(f"Corresponding gaps (N_k values, excluding multiplicities):\n{gaps1}")
    print("The gaps are values from a Poisson distribution, which has unbounded support. Almost surely, the limsup of the gaps is infinity.\n")

    # Set 2
    print("--- Set 2: S = {n^4} ---")
    s2 = get_set_2(20)
    gaps2 = compute_gaps(s2)
    print(f"The first {len(s2)} terms:\n{s2}")
    print(f"Corresponding gaps:\n{gaps2}")
    print("The gaps ( (n+1)^4 - n^4 ) clearly tend to infinity.\n")

    # Set 3
    print("--- Set 3: The set of primes ---")
    s3 = get_set_3(20)
    gaps3 = compute_gaps(s3)
    print(f"The first {len(s3)} primes:\n{s3}")
    print(f"Corresponding gaps:\n{gaps3}")
    print("The gaps between consecutive primes are known to be unbounded.\n")

    # Set 4
    print("--- Set 4: S = {floor((pi/2)^n)} ---")
    s4 = get_set_4(20)
    gaps4 = compute_gaps(s4)
    print(f"The first {len(s4)} terms (unique sorted):\n{s4}")
    print(f"Corresponding gaps:\n{gaps4}")
    print("The gaps tend to infinity as n grows.\n")
    
    print("Conclusion: All four sets have unbounded gaps. Based on Körner's theorem, this is a sufficient condition for the existence of a power series with the desired properties. Therefore, all four sets are valid.")

if __name__ == "__main__":
    main()
