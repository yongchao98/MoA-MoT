import sympy

def find_prime_divisors():
    """
    Finds and prints all prime divisors p for which the number of elements
    of order p in PSL(3, q^2) and PSL(4, q) are equal, where q=12740347.

    The solution is based on finding primes p for which the number of elements
    of that order is zero in both groups. This occurs if ord_p(q) > 4 and
    ord_p(q^2) > 3. The smallest integer d = ord_p(q) that satisfies these
    conditions is d=5.

    Thus, any prime factor of Phi_5(q) = q^4 + q^3 + q^2 + q + 1 is a solution.
    """
    q = 12740347

    # Calculate the value of the 5th cyclotomic polynomial Phi_5(q)
    value_to_factor = q**4 + q**3 + q**2 + q + 1

    # Use sympy's factorint function to find the prime factors
    prime_factors = sympy.factorint(value_to_factor)

    print("The prime divisors p are the factors of q^4 + q^3 + q^2 + q + 1 where q = 12740347.")
    print("Based on the analysis, these primes are:")
    
    # Check if there are any factors found
    if not prime_factors:
        print("No prime factors found (the number might be prime itself).")
        print(value_to_factor)
    else:
        # Print each prime factor
        for p, exponent in prime_factors.items():
            print(p)

# Execute the function
find_prime_divisors()