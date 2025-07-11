import sympy

def solve_problem():
    """
    Finds all prime divisors p of q=12740347 for which the number of elements
    of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347

    print(f"The given number is q = {q}.")

    # Step 1: Find prime divisors of q
    if sympy.isprime(q):
        print(f"The number q = {q} is a prime number.")
        print(f"Therefore, the only prime divisor 'p' of q is p = q itself.")
    else:
        factors = sympy.primefactors(q)
        print(f"The prime divisors of q are: {factors}.")
        print("We will test each of these primes.")
        # This case is not met, but the code structure handles it.

    # Candidate p is q
    p = q
    print(f"\nTesting p = {p}.")

    # Step 2: Explain the counting method for p = characteristic
    print("The prime p is the characteristic of the finite fields F_q and F_{q^2}.")
    print("Elements of order p in PSL(n, k) (where p=char(k)) are the non-identity unipotent elements.")
    print("The number of unipotent elements in GL(n, k) is given by the formula k^(n*(n-1)).")
    print("These all lie in SL(n, k). The number of elements of order p is this value minus 1 (for the identity).")

    # Step 3: Calculate for G1 = PSL(3, q^2)
    n1 = 3
    # The field size k1 is q^2
    print(f"\nFor G1 = PSL(3, q^2):")
    print(f"  n = {n1}, field size k = q^2")
    # Using python's f-string formatting to show the equation
    exponent_n1 = n1 * (n1 - 1)
    final_exponent_n1 = exponent_n1 * 2
    print(f"  Number of elements of order q = (q^2)^({n1}*({n1}-1)) - 1 = (q^2)^{exponent_n1} - 1 = q^{final_exponent_n1} - 1")

    # Step 4: Calculate for G2 = PSL(4, q)
    n2 = 4
    # The field size k2 is q
    print(f"\nFor G2 = PSL(4, q):")
    print(f"  n = {n2}, field size k = q")
    exponent_n2 = n2 * (n2 - 1)
    print(f"  Number of elements of order q = q^({n2}*({n2}-1)) - 1 = q^{exponent_n2} - 1")


    # Step 5: Compare and conclude
    print("\nComparison:")
    print(f"The number of elements of order q in PSL(3, q^2) is q^{final_exponent_n1} - 1.")
    print(f"The number of elements of order q in PSL(4, q) is q^{exponent_n2} - 1.")

    if final_exponent_n1 == exponent_n2:
        print(f"Since {final_exponent_n1} = {exponent_n2}, the number of elements is equal for p = q.")
        print(f"Thus, the only prime divisor p satisfying the condition is p = {q}.")
    else:
        print("The number of elements is not equal.")


solve_problem()