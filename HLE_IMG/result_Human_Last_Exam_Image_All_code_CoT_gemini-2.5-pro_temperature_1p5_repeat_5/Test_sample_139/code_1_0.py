def solve_equation_year():
    """
    This script determines the earliest possible year a given mathematical
    equation could have been written based on the history of its symbols,
    and then rounds that year to the nearest 10.
    """

    # A dictionary mapping key mathematical notations to their year of introduction.
    symbol_introduction_years = {
        "Equals sign '='": 1557,
        "Plus-minus sign '±'": 1631,
        "Superscript exponent 's²'": 1637,
        "Leibniz's differentiation notation 'd/ds'": 1684
    }

    print("Step 1: Identify the introduction year for each key mathematical symbol.")
    for symbol, year in symbol_introduction_years.items():
        print(f"- {symbol} was introduced in {year}.")

    # The first possible year is the latest of these introduction years.
    first_possible_year = max(symbol_introduction_years.values())

    print("\nStep 2: Find the first possible year all symbols were available.")
    # We display the numbers used in the max function.
    years = list(symbol_introduction_years.values())
    print(f"The first possible year is the maximum of the introduction years.")
    print(f"First Possible Year = max({years[0]}, {years[1]}, {years[2]}, {years[3]})")
    print(f"First Possible Year = {first_possible_year}")

    # Step 3: Round the year to the nearest 10.
    base = 10
    # Using the standard rounding method where x.5 rounds up.
    # Python's round() function handles this correctly for rounding to a negative decimal place.
    # For example, round(1684, -1) is 1680.
    rounded_year = int(round(first_possible_year, -1))
    
    print("\nStep 3: Round the year to the nearest 10.")
    print(f"The final calculation is rounding {first_possible_year} to the nearest {base}.")
    # We explicitly show the numbers involved in the rounding calculation.
    print(f"Equation: round({first_possible_year}, -1) = {rounded_year}")
    
    print(f"\nThus, the first possible year this equation could have been written, rounded to the nearest 10 years, is {rounded_year}.")


solve_equation_year()
<<<1680>>>