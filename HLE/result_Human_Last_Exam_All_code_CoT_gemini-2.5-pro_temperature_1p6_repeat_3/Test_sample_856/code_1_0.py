import math

def get_prime_factorization(num):
    """Computes the prime factorization of a number."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n, factors):
    """Calculates Euler's totient function using prime factorization."""
    if n == 1:
        return 1
    result = n
    for p in factors:
        result = result * (p - 1) // p
    return result

def count_solutions_k2_eq_1(factors):
    """
    Counts solutions to k^2 = 1 (mod n) using Chinese Remainder Theorem.
    Number of solutions mod p^a:
    - p=2, a=1: 1
    - p=2, a=2: 2
    - p=2, a>=3: 4
    - p odd: 2
    """
    num_sols = 1
    for p, a in factors.items():
        if p == 2:
            if a == 1:
                num_sols *= 1
            elif a == 2:
                num_sols *= 2
            else: # a >= 3
                num_sols *= 4
        else: # p is an odd prime
            num_sols *= 2
    return num_sols

# The order of the fundamental group
n = math.factorial(10)

# 1. Get the prime factorization of 10!
factors = get_prime_factorization(n)

# 2. Calculate phi(10!)
phi_n = phi(n, factors)

# 3. Calculate the number of solutions to k^2 = 1 (mod 10!)
num_k2_eq_1 = count_solutions_k2_eq_1(factors)

# 4. The number of solutions for k = -k (mod n) and k^2 = -1 (mod n) is 0 for n = 10!
# k=-k -> 2k=0. Requires gcd(k,n)=1, impossible as gcd(n/2, n) != 1 for n > 2.
# k^2=-1. No solutions mod 4, so no solutions mod 10! by CRT.
num_k_eq_neg_k = 0
num_k2_eq_neg_1 = 0

# 5. Apply the formula for the number of non-homeomorphic lens spaces
# The formula is (1/4) * (phi(n) + N(k^2=1) + N(k=-k) + N(k^2=-1))
total_sum = phi_n + num_k2_eq_1 + num_k_eq_neg_k + num_k2_eq_neg_1
num_manifolds = total_sum // 4

# Print the final equation with all numbers
print(f"The number of such manifolds is given by the formula for counting lens spaces L({n}, k):")
print(f"(φ({n}) + N(k²≡1) + N(k≡-k) + N(k²≡-1)) / 4")
print(f"= ({phi_n} + {num_k2_eq_1} + {num_k_eq_neg_k} + {num_k2_eq_neg_1}) / 4")
print(f"= ({phi_n + num_k2_eq_1}) / 4")
print(f"= {num_manifolds}")
