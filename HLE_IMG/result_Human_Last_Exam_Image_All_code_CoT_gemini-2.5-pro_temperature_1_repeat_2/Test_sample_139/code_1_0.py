def solve_equation_year():
    """
    Determines the earliest possible year an equation could have been written
    based on the historical introduction of its symbols.
    """

    # A dictionary of key symbols and their introduction years in a mathematical/typesetting context.
    # The use of '@' as a variable in a formal equation is highly characteristic of modern
    # digital typesetting systems like TeX.
    symbol_years = {
        "Equals sign (=)": 1557,
        "Plus-minus sign (±)": 1631,
        "Newton's dot notation (ṡ)": 1671,
        "At symbol (@) in typeset equation": 1978, # First release of TeX typesetting system
    }

    print("Analyzing the symbols in the equation to find the most recent one:")
    for symbol, year in symbol_years.items():
        print(f"- {symbol}: Year {year}")

    # The earliest possible year is determined by the latest symbol introduction.
    first_possible_year = max(symbol_years.values())

    print(f"\nThe most recently introduced symbol or typesetting practice determines the earliest possible year.")
    print(f"This is the use of the '@' symbol, made practical by systems like TeX, released in {first_possible_year}.")

    # Round the year to the nearest 10 years.
    # Formula: round(year / 10) * 10
    rounded_year = round(first_possible_year / 10) * 10
    
    print(f"\nRounding {first_possible_year} to the nearest 10 years.")
    print(f"The final calculation is: {rounded_year} = round({first_possible_year} / 10) * 10")

solve_equation_year()
<<<1980>>>