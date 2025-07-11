def solve_equation_year():
    """
    Calculates the first possible year an equation could have been written
    based on the history of its symbols.
    """
    # Historical dates for the introduction of key symbols into mathematical/typesetting use.
    # The use of '@' as a variable is tied to its availability on standard keyboards.
    symbol_dates = {
        'equals_sign': 1557,
        'plus_minus_sign': 1631,
        'exponent_notation': 1637,
        'differential_notation': 1684,
        'at_symbol_on_keyboard': 1874  # Sholes & Glidden typewriter
    }

    # The first possible year is determined by the most recent symbol invention.
    first_possible_year = max(symbol_dates.values())

    # The problem asks to round the year to the nearest 10 years.
    # The round() function with a negative number of digits rounds to the nearest 10, 100, etc.
    # round(1874, -1) correctly rounds to 1870.
    rounded_year = int(round(first_possible_year, -1))

    print(f"The analysis of symbols suggests the following timeline:")
    print(f"- Plus-minus sign ('±'): {symbol_dates['plus_minus_sign']}")
    print(f"- Differential notation ('d/ds'): {symbol_dates['differential_notation']}")
    print(f"- '@' symbol on a standard keyboard: {symbol_dates['at_symbol_on_keyboard']}")
    print(f"\nThe limiting factor is the '@' symbol's availability on a typewriter, making the first possible year {first_possible_year}.")
    print(f"Rounding {first_possible_year} to the nearest 10 years gives {rounded_year}.")

solve_equation_year()
<<<1870>>>