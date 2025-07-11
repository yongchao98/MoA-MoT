import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function."""
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

def solve_asymptotic_exponents():
    """
    Calculates the exponents alpha and beta in the asymptotic formula for |A(X)|.
    """
    # The order of the character group is K=12.
    K = 12

    # The exponent beta is derived from a value 'w' related to the structure
    # of the character group (Z/12Z)*.
    # w = (1/phi(K)) * sum_{a in (Z/KZ)*} gcd(a-1, K)

    phi_K = phi(K)
    
    # Find the elements of the multiplicative group of integers modulo K, (Z/KZ)*.
    group_elements = [a for a in range(1, K) if gcd(a, K) == 1]

    # Calculate the sum of gcd(a-1, K) for each element a in the group.
    sum_of_gcds = 0
    print("The problem is about characters of order dividing 12.")
    print(f"We analyze the group (Z/12Z)*, which has phi(12) = {phi_K} elements.")
    print(f"The elements are: {group_elements}")
    print("\nWe compute gcd(a-1, 12) for each element 'a':")
    for a in group_elements:
        current_gcd = gcd(a - 1, K)
        print(f"  a = {a:2d}, gcd({a-1:2d}, 12) = {current_gcd}")
        sum_of_gcds += current_gcd

    # Calculate w.
    w = sum_of_gcds / phi_K
    print(f"\nThe sum of these gcds is {sum_of_gcds}.")
    print(f"The value w = (sum of gcds) / phi(12) = {sum_of_gcds} / {phi_K} = {w:.1f}")

    # The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta.
    # From analytic number theory, alpha is 1 and beta is w - 2.
    alpha = 1
    beta = w - 2

    print("\nThe exponents in the asymptotic formula are derived from w:")
    print(f"alpha = {int(alpha)}")
    print(f"beta = w - 2 = {int(w)} - 2 = {int(beta)}")

    sum_alpha_beta = alpha + beta
    print("\nThe final equation is |A(X)| ~ c * X^1 * log(X)^4")
    print(f"The sum of the integers alpha and beta is:")
    print(f"alpha + beta = {int(alpha)} + {int(beta)} = {int(sum_alpha_beta)}")

solve_asymptotic_exponents()