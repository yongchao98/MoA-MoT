def solve():
    """
    Solves the problem by calculating a(n) mod p for two given primes p.
    a(n) is the number of domino tilings of a 3x(2n) rectangle.
    """

    # Function to calculate a(n) using the recurrence a(n) = 4*a(n-1) - a(n-2)
    def calculate_a(n):
        if n == 0:
            return 1
        a0, a1 = 1, 3
        if n == 1:
            return a1
        
        # Iteratively calculate a(n)
        for i in range(2, n + 1):
            a_next = 4 * a1 - a0
            a0, a1 = a1, a_next
        return a1

    results = []

    # Case 1: p = 50051
    p1 = 50051
    # The argument to a(n) simplifies to 5
    n1_effective = 5
    # The expression to be calculated
    equation_str1 = f"a({p1}^4 + 4*{p1}^3 - 5*{p1}^2 - 3*{p1} + 8)"
    
    # Calculate a(5)
    result1 = calculate_a(n1_effective)
    results.append(result1)
    
    print(f"For p = {p1}:")
    print(f"{equation_str1} mod {p1} simplifies to a({n1_effective}) mod {p1}.")
    a_values_str = []
    a_minus_2, a_minus_1 = 1, 3
    a_values_str.append("a(0) = 1")
    a_values_str.append("a(1) = 3")
    for i in range(2, n1_effective + 1):
        a_current = 4 * a_minus_1 - a_minus_2
        a_values_str.append(f"a({i}) = 4*{a_minus_1} - {a_minus_2} = {a_current}")
        a_minus_2, a_minus_1 = a_minus_1, a_current
    print("Calculation: " + ", ".join(a_values_str))
    print(f"Result for p={p1} is {result1}\n")


    # Case 2: p = 50069
    p2 = 50069
    # The argument to a(n) simplifies to 8, with a negative sign
    n2_effective = 8
    # The expression to be calculated
    equation_str2 = f"a({p2}^4 + 4*{p2}^3 - 5*{p2}^2 - 3*{p2} + 8)"
    
    # Calculate a(8)
    a_8 = calculate_a(n2_effective)
    result2 = (-a_8) % p2
    results.append(result2)

    print(f"For p = {p2}:")
    print(f"{equation_str2} mod {p2} simplifies to -a({n2_effective}) mod {p2}.")
    print(f"First, we calculate a(8) = {a_8}.")
    print(f"Then, the result is (-{a_8}) mod {p2} = {result2}.")
    print(f"Result for p={p2} is {result2}\n")

    # Final combined answer
    print("Final comma-separated answers:")
    print(f"{results[0]},{results[1]}")

solve()