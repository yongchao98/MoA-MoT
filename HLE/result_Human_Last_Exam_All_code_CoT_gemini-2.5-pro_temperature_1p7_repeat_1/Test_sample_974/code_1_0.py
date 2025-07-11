import sympy

def find_primes():
    """
    Solves the user's question by finding all prime divisors p of q=12740347
    for which the number of elements of order p in both PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347

    # Step 1: Find the prime divisors of q.
    # We use the sympy library to get the prime factors.
    prime_divisors = sympy.primefactors(q)

    # Step 2: Since q is prime, the only prime divisor is q itself.
    # We must check if the condition holds for p = q.
    # The number of elements of order p=char(k) in PSL(n, k) is k^(n(n-1)) - 1.

    # Step 3: Calculate the number of elements for G1 = PSL(3, q^2)
    # Here n=3, k=q^2.
    n1 = 3
    # The exponent is n1 * (n1 - 1) = 3 * 2 = 6. The base is q^2.
    # Total exponent for q is 2 * 6 = 12.
    final_exponent1 = 12

    # Step 4: Calculate the number of elements for G2 = PSL(4, q)
    # Here n=4, k=q.
    n2 = 4
    # The exponent is n2 * (n2 - 1) = 4 * 3 = 12. The base is q.
    # Total exponent for q is 12.
    final_exponent2 = 12

    # Step 5: Explain the findings and present the result.
    print(f"Let q = {q}.")
    print("The task is to find all prime divisors p of q such that:")
    print("Number of elements of order p in PSL(3, q^2) = Number of elements of order p in PSL(4, q)")
    print("-" * 30)

    print(f"First, we find the prime factorization of q = {q}.")
    if len(prime_divisors) == 1 and prime_divisors[0] == q:
        print(f"The number q = {q} is a prime number.")
        print(f"Therefore, the only prime divisor p of q to check is p = {q}.")
    else:
        print(f"The prime divisors of q are: {prime_divisors}")

    print("\nSince p=q is the characteristic of the fields, we count unipotent elements.")
    print("The number of elements of order p in PSL(n, k) (where char(k)=p) is k^(n*(n-1)) - 1.")

    print("\nFor G1 = PSL(3, q^2):")
    print(f"The number of elements of order q is (q^2)^({n1}*({n1}-1)) - 1 = q^{final_exponent1} - 1.")

    print("\nFor G2 = PSL(4, q):")
    print(f"The number of elements of order q is q^({n2}*({n2}-1)) - 1 = q^{final_exponent2} - 1.")

    print("\nNow, we form the final equation by equating the two quantities:")
    print(f"{q}^{final_exponent1} - 1 = {q}^{final_exponent2} - 1")
    
    print("\nThis equality is always true, so the condition holds for p = q.")
    print("The list of all such prime divisors p is therefore:")
    for p in prime_divisors:
        print(p)

find_primes()