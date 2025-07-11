def solve_piazza_mystery():
    """
    This function identifies and prints the architectural feature
    removed from Piazza della Rotonda in 1950.
    """
    # The year the feature was removed.
    year = 1950
    
    # The feature that was removed.
    removed_feature = "tram tracks"

    # Print the answer to the user's question.
    print(f"The unique feature in the Piazza della Rotonda that was removed in {year} was the set of {removed_feature}.")
    print("The tram line, which ran through the square, was considered disruptive to the historic site and was removed.")

    # As requested, here are the numbers from the final 'equation' (the year of removal).
    print(f"\nThe numbers from the year {year} are:")
    year_str = str(year)
    for digit in year_str:
        print(digit)

# Execute the function
solve_piazza_mystery()