def solve_piazza_mystery():
    """
    Calculates the year an architectural feature was removed from
    Piazza della Rotonda and identifies the feature.
    """
    # Define the numerical components to form the year 1950
    num1 = 1000
    num2 = 900
    num3 = 50

    # Calculate the year
    year_removed = num1 + num2 + num3

    # As requested, output the full equation with each number
    print("Calculating the significant year for the Piazza della Rotonda:")
    print(f"{num1} + {num2} + {num3} = {year_removed}")

    # Print the answer to the user's question
    print(f"\nIn the year {year_removed}, the feature that was removed from the Piazza della Rotonda was its surrounding iron railings.")

solve_piazza_mystery()