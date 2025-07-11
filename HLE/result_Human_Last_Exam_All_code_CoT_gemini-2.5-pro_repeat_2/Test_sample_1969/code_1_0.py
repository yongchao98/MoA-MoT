def calculate_and_show_shapley_value():
    """
    Asks the user for n and k, then calculates and displays the Shapley value
    (fair share) for player p_k based on the derived formula.
    """
    try:
        n_str = input("Enter the total number of people (n > 1): ")
        n = int(n_str)
        if n <= 1:
            print("Error: n must be an integer greater than 1.")
            return

        k_str = input(f"Enter the player index k (1 <= k <= {n}): ")
        k = int(k_str)
        if not (1 <= k <= n):
            print(f"Error: k must be an integer between 1 and {n}.")
            return
    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # Calculate S1 = sum of first n integers
    s1 = n * (n + 1) // 2
    
    # Calculate S2 = sum of first n squares
    s2 = n * (n + 1) * (2 * n + 1) // 6
    
    # Pre-calculate powers of S1 for use in the formula
    s1_squared = s1 ** 2
    s1_cubed = s1 ** 3

    # The simplified formula for the Shapley value c_k is:
    # c_k = k*S1*S2 + k*S1^3 - k^2*S1^2
    
    # Calculate each term of the formula
    term1 = k * s1 * s2
    term2 = k * s1_cubed
    term3 = (k ** 2) * s1_squared
    
    # Calculate the final value for c_k
    c_k = term1 + term2 - term3
    
    # Print the equation with the calculated numbers
    print("\n-------------------------------------------")
    print(f"Calculating the fair share for player p_{k} when n = {n}:")
    print("-------------------------------------------")
    print(f"The general formula for c_k is: c_k = k*S1*S2 + k*S1^3 - k^2*S1^2")
    print("\nFirst, we calculate the values of S1 and S2:")
    print(f"S1 = Sum(1..{n}) = {s1}")
    print(f"S2 = Sum(1^2..{n}^2) = {s2}")
    
    print("\nNext, we substitute these values into the formula:")
    print(f"c_{k} = {k} * ({s1}) * ({s2}) + {k} * ({s1})^3 - ({k})^2 * ({s1})^2")
    print(f"c_{k} = {k} * {s1} * {s2} + {k} * {s1_cubed} - {k**2} * {s1_squared}")
    
    print("\nCalculating each term:")
    print(f"Term 1: {k} * {s1} * {s2} = {term1}")
    print(f"Term 2: {k} * {s1_cubed} = {term2}")
    print(f"Term 3: {k**2} * {s1_squared} = {term3}")
    
    print("\nFinally, combining the terms:")
    print(f"c_{k} = {term1} + {term2} - {term3}")
    print(f"\nThe exact amount of money for player p_{k} is: ${c_k}")
    print("-------------------------------------------\n")

# Execute the function
calculate_and_show_shapley_value()