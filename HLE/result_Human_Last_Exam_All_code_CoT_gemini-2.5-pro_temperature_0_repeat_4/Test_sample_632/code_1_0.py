def solve_piazza_history():
    """
    This function provides the answer to the historical question about the
    Piazza della Rotonda and displays the relevant numbers as requested.
    """
    # The unique feature that was removed from the square
    feature = "The tramway rails"

    # The year the feature was removed by the city administration
    year_of_removal = 1950

    print(f"The unique architectural feature removed from the Piazza della Rotonda in Rome was: {feature}.")
    print(f"This was done by the city administration in the year {year_of_removal}.")
    print("\nAs requested, here are the individual numbers from the final 'equation' (the year):")

    # Output each number (digit) from the year
    for digit in str(year_of_removal):
        print(digit)

# Execute the function to provide the answer
solve_piazza_history()