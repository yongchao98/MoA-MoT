def solve():
    """
    This function calculates the number of starting positions in the described game
    where the bitwise XOR value of the piles' Grundy values is one or two.
    """
    # Per the problem statement, n > 200 and t > 0.
    # We will use example values for demonstration.
    n = 250
    t = 5

    print(f"Solving for n = {n} piles and parameter t = {t}.")
    print("-" * 50)
    print("This is an impartial game, so we use the Sprague-Grundy theorem.")
    print("The game's outcome depends on the nim-sum (XOR sum) of the Grundy values of the piles.")
    print("Let G_a be the Grundy value for a pile of 2^a stones.")
    print("-" * 50)

    print("\nStep 1: Determining the Grundy values G_a = g(2^a)")
    G = {}
    G[0] = 1  # G_0 = g(2^0) = g(1) = mex{g(0)} = 1
    G[1] = 0  # G_1 = g(2^1) = g(2) = mex{g(1)} = mex{1} = 0
    G[2] = 2  # G_2 = g(2^2) = g(4) = mex{g(2), g(1)} = mex{0, 1} = 2
    print("G_0 = 1")
    print("G_1 = 0")
    print("G_2 = 2")
    print("For a > 2, the recurrence is G_a = mex{G_{a-1}, G_{a-2}, G_{a-3}}.")
    print("This results in a period-4 sequence for a >= 3: (3, 1, 0, 2, ...)")

    print("-" * 50)
    print(f"\nStep 2: Counting choices for exponents a_i in [1, {4*t+2}] for each Grundy value")
    N0 = t + 1
    N1 = t
    N2 = t + 1
    N3 = t
    total_choices_per_pile = 4 * t + 2
    
    print(f"Number of choices for a_i yielding G_a=0 is N0 = t + 1 = {N0}")
    print(f"Number of choices for a_i yielding G_a=1 is N1 = t = {N1}")
    print(f"Number of choices for a_i yielding G_a=2 is N2 = t + 1 = {N2}")
    print(f"Number of choices for a_i yielding G_a=3 is N3 = t = {N3}")
    print(f"Total number of choices for each a_i is N0+N1+N2+N3 = {total_choices_per_pile}.")
    
    print("-" * 50)
    print("\nStep 3: Calculating the number of starting positions with nim-sum 1 or 2.")
    print("Let Num(S) be the number of positions with nim-sum S.")
    print(f"The number of positions with nim-sum=1 is given by the formula: (({total_choices_per_pile})^n - 2^n) / 4")
    print(f"The number of positions with nim-sum=2 is given by the formula: (({total_choices_per_pile})^n + 2^n) / 4")
    
    base = total_choices_per_pile
    term1 = pow(base, n)
    term2 = pow(2, n)

    num1 = (term1 - term2) // 4
    num2 = (term1 + term2) // 4
    total = num1 + num2
    
    print("\nCalculating for nim-sum = 1:")
    print(f"({base}^{n} - 2^{n}) / 4")
    print(f"= ({term1} - {term2}) / 4")
    print(f"= {num1}")

    print("\nCalculating for nim-sum = 2:")
    print(f"({base}^{n} + 2^{n}) / 4")
    print(f"= ({term1} + {term2}) / 4")
    print(f"= {num2}")

    print("\nTotal number of starting positions with nim-sum 1 or 2 is the sum:")
    print(f"{num1} + {num2} = {total}")

    print("-" * 50)
    print("\nFinal Result Verification")
    print("The total number simplifies to the final formula: (4*t + 2)^n / 2")
    final_result = pow(4*t+2, n) // 2
    print(f"Final Calculation: ({base})^{n} / 2 = {term1} / 2 = {final_result}")

if __name__ == '__main__':
    solve()