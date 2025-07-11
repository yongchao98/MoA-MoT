import math

def get_prime_divisors(n):
    """
    Returns a sorted list of unique prime divisors of a number n.
    """
    factors = set()
    d = 2
    temp_n = n
    # Check for factors of 2
    if temp_n % d == 0:
        factors.add(d)
        while temp_n % d == 0:
            temp_n //= d
    # Check for odd factors
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

def solve_and_print():
    """
    Solves the problem and prints the reasoning and the result.
    """
    q = 12740347

    # Step 1: Find the prime divisors of q. The problem is interpreted as
    # finding prime divisors p OF q.
    prime_divs_of_q = get_prime_divisors(q)

    # Step 2: For each prime divisor p (which must be q), check the condition.
    # A prime divisor p of q must be q itself, since q is prime.
    # For this prime p=q, it is the characteristic of the fields F_q and F_{q^2}.
    # The number of elements of order p in PSL(n, r) where p = char(F_r) is
    # given by N_p = r^(n*(n-1)) - 1.

    result_primes = []
    for p in prime_divs_of_q:
        # Check for G1 = PSL(3, q^2), where n=3, r=q^2
        n1 = 3
        exp1_base = n1 * (n1 - 1)
        exp1_final = exp1_base * 2 # from (q^2)^(n*(n-1))

        # Check for G2 = PSL(4, q), where n=4, r=q
        n2 = 4
        exp2_final = n2 * (n2 - 1)

        # If the final exponents are equal, the numbers of elements are equal.
        if exp1_final == exp2_final:
            result_primes.append(p)

    # Step 3: Print the explanation and the final answer.
    print(f"Let q = {q}.")
    print("We need to find prime divisors p of q for which the number of elements of order p")
    print("in PSL(3, q^2) and PSL(4, q) are equal.")
    print(f"\nThe prime divisors of q={q} are: {prime_divs_of_q}")
    print("-" * 50)

    if not result_primes:
        print("No such prime divisors found.")
    else:
        for p in result_primes:
            print(f"Checking for prime divisor p = {p}:")
            print("Since p is a divisor of q (a prime), p must be equal to q.")
            print("p is the characteristic of the fields F_q and F_{q^2}.")
            print("The number of elements of order p in PSL(n, r) (where p is the characteristic)")
            print("is given by the formula: r^(n*(n-1)) - 1.\n")

            # Details for PSL(3, q^2)
            n1, r_str1 = 3, f"{p}^2"
            exp1_base = n1 * (n1 - 1)
            exp1_final = exp1_base * 2
            print(f"For PSL({n1}, {p}^2):")
            print(f"n = {n1}, r = {r_str1}")
            print(f"Number of elements of order {p} = ({r_str1})^({n1}*({n1}-1)) - 1 = ({r_str1})^{exp1_base} - 1 = {p}^{exp1_final} - 1")

            # Details for PSL(4, q)
            n2, r_str2 = 4, f"{p}"
            exp2_final = n2 * (n2 - 1)
            print(f"\nFor PSL({n2}, {p}):")
            print(f"n = {n2}, r = {r_str2}")
            print(f"Number of elements of order {p} = {r_str2}^({n2}*({n2}-1)) - 1 = {p}^{exp2_final} - 1")

            print("\nComparing the two final equations:")
            print(f"Equation 1: N_p(PSL(3, q^2)) = {p}^{exp1_final} - 1")
            print(f"Equation 2: N_p(PSL(4, q)) = {p}^{exp2_final} - 1")
            if exp1_final == exp2_final:
                print(f"Since {exp1_final} = {exp2_final}, the numbers of elements are equal. So, p={p} is a solution.")

        print("-" * 50)
        print("The list of all such prime divisors is:")
        for p in result_primes:
            print(p)

solve_and_print()