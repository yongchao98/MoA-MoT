import math

def get_prime_divisors(n):
    """
    This function returns a list of unique prime divisors of a given number n.
    """
    factors = set()
    d = 2
    temp_n = n
    # Check for divisibility by 2
    if temp_n % 2 == 0:
        factors.add(2)
        while temp_n % 2 == 0:
            temp_n //= 2
    # Check for divisibility by odd numbers
    d = 3
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 2
    if temp_n > 1:
        factors.add(temp_n)
    return sorted(list(factors))

def solve_problem():
    """
    Solves the problem by finding prime divisors p of q for which the number of
    elements of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347
    print(f"The given number is q = {q}")

    prime_divisors_of_q = get_prime_divisors(q)
    print(f"The prime divisors of q are: {prime_divisors_of_q}")
    print("-" * 30)

    solution_primes = []

    # For each prime divisor p, check the condition.
    # The formula used is valid when p is the characteristic of the field.
    # The fields are F_q and F_{q^2}, and their characteristic is q.
    # Since q is prime, its only prime divisor is p=q, which is the characteristic.
    for p in prime_divisors_of_q:
        if p == q:
            # Case 1: PSL(3, q^2)
            # The number of elements of order p=q is f^(n(n-1)) - 1,
            # where n=3 and f = q^2.
            n1 = 3
            # The formula gives: (q^2)^(3 * (3-1)) - 1 = q^12 - 1

            # Case 2: PSL(4, q)
            # The number of elements of order p=q is f^(n(n-1)) - 1,
            # where n=4 and f = q.
            n2 = 4
            # The formula gives: q^(4 * (4-1)) - 1 = q^12 - 1

            # We check if the exponents are equal.
            exponent1 = 2 * n1 * (n1 - 1)
            exponent2 = n2 * (n2 - 1)

            if exponent1 == exponent2:
                solution_primes.append(p)
                print(f"Checking for prime divisor p = {p}...")
                print("This prime is the characteristic of the fields F_q and F_{q^2}.")
                
                print("\nFor PSL(3, q^2):")
                print(f"The number of elements of order {p} is (q^2)^(n*(n-1)) - 1 with n=3, q={q}")
                print(f"The equation is: ({q}^2)^({n1}*({n1}-1)) - 1 = {q}^({2*n1*(n1-1)}) - 1 = {q}^{exponent1} - 1")

                print("\nFor PSL(4, q):")
                print(f"The number of elements of order {p} is q^(n*(n-1)) - 1 with n=4, q={q}")
                print(f"The equation is: {q}^({n2}*({n2}-1)) - 1 = {q}^{exponent2} - 1")

                print(f"\nComparing the two results: {q}^{exponent1} - 1 == {q}^{exponent2} - 1.")
                print("The numbers of elements are equal.")
            else:
                print(f"For p = {p}, the numbers of elements are not equal.")

    print("-" * 30)
    print("The prime divisor(s) p for which the condition holds is/are:")
    if not solution_primes:
        print("None")
    else:
        for p in solution_primes:
            print(p)

solve_problem()