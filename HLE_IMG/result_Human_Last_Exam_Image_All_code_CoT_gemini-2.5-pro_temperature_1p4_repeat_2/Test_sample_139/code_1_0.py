import math

def solve_year_puzzle():
    """
    Determines the first possible year an equation could have been written
    based on the historical introduction of its symbols.
    """
    # Introduction year of the Plus-Minus symbol (±)
    year_plus_minus = 1631

    # Year by which calculus notations (both Newton's and Leibniz's) were established.
    # Newton's Method of Fluxions was published in 1736.
    year_calculus = 1736

    # The year the '@' symbol became widely available on a standard keyboard.
    # The Underwood No. 5 typewriter, introduced in 1900, made it common.
    # This is the limiting factor for the symbol's use in a typeset paper.
    year_at_sign = 1900

    # The first possible year is the latest of these dates, as all symbols
    # must have been available.
    years = [year_plus_minus, year_calculus, year_at_sign]
    first_possible_year = max(years)

    # Round the year to the nearest 10 years.
    # A negative second argument to round() rounds to powers of 10.
    rounded_year = int(round(first_possible_year, -1))

    # Output the explanation and the final equation for clarity.
    print("To find the earliest possible year, we must find the latest introduction date among the symbols used.")
    print(f"The '±' symbol was introduced around: {year_plus_minus}")
    print(f"Calculus notation was established by: {year_calculus}")
    print(f"The '@' symbol became available on standard typewriters around: {year_at_sign}")
    print("\nThe final calculation to find the limiting year is:")
    print(f"max({year_plus_minus}, {year_calculus}, {year_at_sign}) = {first_possible_year}")

    print(f"\nRounding {first_possible_year} to the nearest 10 years gives the final answer.")
    print(f"Answer: {rounded_year}")

solve_year_puzzle()
<<<1900>>>