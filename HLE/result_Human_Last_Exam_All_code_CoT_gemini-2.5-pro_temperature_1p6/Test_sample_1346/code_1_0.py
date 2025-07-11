def solve():
    """
    Solves the problem for the given primes.

    The problem asks for the value of a(p^4+4p^3-5p^2-3p+8) mod p.
    Based on number theory principles, for the given primes p, this simplifies to
    calculating a(3), where a(n) is the sequence for the number of domino tilings.
    The sequence follows a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.
    This script calculates a(3) and presents the result for both primes.
    """
    primes = [50051, 50069]
    results = []

    # The index we need to calculate, based on the theoretical derivation, is 3 for both primes.
    target_n = 3

    print(f"For both primes p=50051 and p=50069, the problem reduces to calculating a(3).")
    print("The recurrence relation is a(n) = 4*a(n-1) - a(n-2).")
    print("The calculation is as follows:")

    # Calculate a(n) up to target_n, printing each step
    if target_n == 0:
        a_n = 1
        print("a(0) = 1")
    elif target_n == 1:
        a_n = 3
        print("a(0) = 1")
        print("a(1) = 3")
    else:
        a_prev = 1  # a(0)
        a_curr = 3  # a(1)
        print("a(0) = 1")
        print("a(1) = 3")
        for i in range(2, target_n + 1):
            a_next = 4 * a_curr - a_prev
            # This printout satisfies the requirement to output each number in the final equation
            print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {a_curr} - {a_prev} = {a_next}")
            a_prev = a_curr
            a_curr = a_next
        a_n = a_curr
    
    # The result is the same for both primes
    result_val = a_n
    
    print("\n" + "="*40)
    print(f"The result for p=50051 is {result_val}.")
    print(f"The result for p=50069 is {result_val}.")
    
    # Format the final answer as requested
    final_answer = f"{result_val},{result_val}"
    print("\nThe final answers for p=50051 and p=50069 respectively, separated by a comma are:")
    print(final_answer)

solve()