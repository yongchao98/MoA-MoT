def solve_game_positions():
    """
    Calculates the number of starting positions for a game on piles of stones
    where the bitwise XOR sum of the Grundy values of the piles is 1 or 2.

    The user is prompted to enter the number of piles (n) and a parameter (t),
    where n > 200 and t > 0.
    """
    try:
        # Prompt the user to enter the values for n and t
        n_str = input("Enter the number of piles n (n > 200): ")
        t_str = input("Enter the integer t (t > 0): ")

        n = int(n_str)
        t = int(t_str)

        if n <= 200 or t <= 0:
            print("Invalid input: n must be greater than 200 and t must be greater than 0.")
            return

        # The number of starting positions is given by the formula:
        # 2^(n-1) * ((2t+1)^n - (-1)^n)
        
        # Calculate each part of the formula. Python's integers handle large numbers.
        term1 = pow(2, n - 1)
        
        base_term2 = 2 * t + 1
        term2 = pow(base_term2, n)
        
        term3 = pow(-1, n)
        
        # Final calculation
        result = term1 * (term2 - term3)

        # Print the breakdown of the calculation as requested
        print("\nThe number of starting positions with a total Grundy value of 1 or 2 can be calculated as follows:")
        print(f"n = {n}")
        print(f"t = {t}")
        print("\nFormula: 2^(n-1) * ((2t+1)^n - (-1)^n)")
        
        # Displaying the intermediate computed values
        print("\nBreakdown of the calculation:")
        print(f"2^(n-1) = 2^({n}-1) = {term1}")
        print(f"(2t+1)^n = (2*{t}+1)^{n} = {base_term2}^{n} = {term2}")
        print(f"(-1)^n = (-1)^{n} = {term3}")
        
        # Displaying the final result in equation form
        print(f"\nFinal calculation: {term1} * ({term2} - ({term3}))")
        print(f"Result: {result}")

    except ValueError:
        print("Invalid input. Please enter valid integers for n and t.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the function
solve_game_positions()
