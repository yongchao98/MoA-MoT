import math

def is_prime(n):
    """
    Checks if a number n is prime.
    """
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_primes_for_equal_orders():
    """
    Solves the problem for the given q.
    """
    q = 12740347

    print(f"Given q = {q}.")
    print("The task is to find all prime divisors p of q such that the number of elements of order p is the same in PSL(3,q^2) and PSL(4,q).")
    
    print("\nStep 1: Find prime divisors of q.")
    if is_prime(q):
        prime_divisors = [q]
        print(f"q = {q} is a prime number.")
        print(f"The only prime divisor to check is p = {q}.")
    else:
        # This part is not needed for the given q, but would be for a composite q.
        print(f"q = {q} is composite. Finding its prime divisors...")
        # Factorization logic would go here.
        return

    # The only prime p to check is q itself.
    p = q
    
    print(f"\nStep 2: Compare the number of elements of order p = q in both groups.")
    print("For elements of order q (the characteristic of the field), the number in PSL(n, r) (with r a power of q) is given by r^(n*(n-1)) - 1.")

    # Group 1: PSL(3, q^2)
    n1 = 3
    # The field is F_q^2, so the base 'r' in the formula is q^2.
    r_base_str1 = "q^2"
    exp1_n = n1 * (n1 - 1)
    exp1_final = 2 * exp1_n
    
    # Group 2: PSL(4, q)
    n2 = 4
    # The field is F_q, so the base 'r' is q.
    r_base_str2 = "q"
    exp2_final = n2 * (n2 - 1)
    
    print("\nFor G1 = PSL(3, q^2):")
    print(f"  n = {n1}, r = q^2")
    print(f"  Number of elements = (q^2)^({n1}*({n1}-1)) - 1 = (q^2)^{exp1_n} - 1 = q^{exp1_final} - 1")
    
    print("\nFor G2 = PSL(4, q):")
    print(f"  n = {n2}, r = q")
    print(f"  Number of elements = q^({n2}*({n2}-1)) - 1 = q^{exp2_final} - 1")

    print("\nStep 3: The Final Equation.")
    print(f"We need to check if q^{exp1_final} - 1 = q^{exp2_final} - 1.")
    print(f"This is true if and only if the exponents are equal: {exp1_final} = {exp2_final}.")
    
    if exp1_final == exp2_final:
        print(f"The exponents are indeed equal.")
        print(f"\nConclusion: The condition holds for the prime divisor p = {p}.")
        print("The list of all such primes is:")
        print(p)
        return [p]
    else:
        print(f"The exponents are not equal.")
        print("\nConclusion: The condition does not hold for p = {p}.")
        print("There are no such primes.")
        return []

find_primes_for_equal_orders()