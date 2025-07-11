def solve_prime_divisors():
    """
    This script finds all prime divisors p of q=12740347 for which the number
    of elements of order p in both PSL(3,q^2) and PSL(4,q) are equal.
    """
    q = 12740347

    def get_prime_divisors(n):
        """Returns a list of prime divisors of a number n."""
        factors = set()
        d = 2
        temp_n = n
        while d * d <= temp_n:
            if temp_n % d == 0:
                factors.add(d)
                while temp_n % d == 0:
                    temp_n //= d
            d += 1
        if temp_n > 1:
            factors.add(temp_n)
        return sorted(list(factors))

    prime_divisors = get_prime_divisors(q)

    print(f"The number is q = {q}.")
    print(f"The prime divisors of q are: {prime_divisors}.")
    print("-" * 20)

    solution_primes = []

    for p in prime_divisors:
        print(f"Checking for prime p = {p}:")
        
        # Analysis for PSL(3, q^2)
        # For n=3, max Jordan block size is 3.
        # If p >= 3, all non-identity unipotents have order p.
        num_g1_expr = "q^(12) - 1"
        
        # Analysis for PSL(4, q)
        # For n=4, max Jordan block size is 4.
        # If p >= 5, all non-identity unipotents have order p.
        if p >= 5:
            num_g2_expr = "q^(12) - 1"
            are_equal = True
        else: # p must be 3, since 2 is not a divisor of q
            num_g2_expr = "q^(12) - 1 - |C_(4)|"
            are_equal = False
        
        print(f"Number of elements of order {p} in PSL(3, q^2) = {num_g1_expr}")
        print(f"Number of elements of order {p} in PSL(4, q) = {num_g2_expr}")

        if are_equal:
            print("The numbers are equal.")
            solution_primes.append(p)
        else:
            print("The numbers are not equal, because for PSL(4,q) we must exclude unipotent classes with Jordan blocks of size > 3.")
        print("-" * 20)

    print("\nFinal Answer:")
    print("The prime divisors p for which the number of elements of order p are equal in both groups are:")
    print(solution_primes)

solve_prime_divisors()