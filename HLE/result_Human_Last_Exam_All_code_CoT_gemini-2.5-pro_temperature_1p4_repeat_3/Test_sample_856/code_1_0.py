import math

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary.
    """
    factors = {}
    d = 2
    temp = n
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def euler_totient(factors):
    """
    Calculates phi(n) from its prime factorization.
    phi(n) = n * product_{p|n} (1 - 1/p)
    """
    n = 1
    for p, a in factors.items():
        n *= (p**a)
    
    phi = n
    for p in factors.keys():
        phi = phi // p * (p - 1)
    return phi

def count_s2_solutions(factors):
    """
    Counts the number of solutions to x^2 = 1 (mod n) from prime factorization of n.
    """
    s2 = 1
    for p, a in factors.items():
        if p == 2:
            if a == 1:
                s2 *= 1
            elif a == 2:
                s2 *= 2
            else: # a >= 3
                s2 *= 4
        else: # odd prime
            s2 *= 2
    return s2

# Calculate for n = 10!
n = math.factorial(10)

# 1. Get prime factorization of 10!
# 10! = 10*9*8*7*6*5*4*3*2*1 = (2*5)*(3^2)*(2^3)*7*(2*3)*5*(2^2)*3*2 = 2^8 * 3^4 * 5^2 * 7^1
factors = {2: 8, 3: 4, 5: 2, 7: 1}

# 2. Calculate phi(10!)
phi_n = euler_totient(factors)

# 3. Calculate S_2(10!)
s2_n = count_s2_solutions(factors)

# 4. Calculate the number of manifolds
# We assume we only count lens spaces (for cyclic fundamental group)
count = (phi_n + s2_n) // 2

# Print the final equation as requested
print(f"{count} = ({phi_n} + {s2_n}) / 2")