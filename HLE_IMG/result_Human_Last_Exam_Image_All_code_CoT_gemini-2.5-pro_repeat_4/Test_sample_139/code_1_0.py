def solve_equation_year():
    """
    This function determines the first possible year an equation could have been
    written based on the introduction of its symbols, and rounds it to the
    nearest 10 years.
    """
    # A dictionary mapping symbols to their year of introduction or popularization.
    symbol_introduction_years = {
        '=': 1557,
        '±': 1631,
        '²': 1637,
        'd/ds': 1675,
        '@': 1874  # Year the '@' symbol was included on the first commercial typewriter.
    }

    # The first possible year is the latest year any symbol was introduced.
    first_possible_year = max(symbol_introduction_years.values())

    # Round the year to the nearest 10.
    # We do this by dividing by 10, rounding to the nearest whole number,
    # and then multiplying by 10.
    final_year = round(first_possible_year / 10) * 10

    print("The most recently introduced symbol in the equation is the '@' symbol.")
    print(f"It became widely available for typesetting with the first commercial typewriter in {symbol_introduction_years['@']}.")
    print(f"Therefore, the earliest possible year the equation could have been written is {first_possible_year}.")
    print(f"Rounding this to the nearest 10 years gives the final answer.")
    print(f"Equation: round({first_possible_year} / 10) * 10 = {final_year}")


solve_equation_year()
<<<1870>>>