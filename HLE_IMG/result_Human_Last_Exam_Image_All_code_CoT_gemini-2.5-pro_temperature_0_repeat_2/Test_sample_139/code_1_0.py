import math

def solve_equation_year():
    """
    Calculates the earliest possible year an equation could have been written
    based on the introduction dates of its symbols.
    """
    # Historical introduction years for each symbol in the equation.
    # The year for '@' is based on its inclusion on commercial typewriters.
    symbol_years = {
        "Equals sign '='": 1557,
        "Plus-minus sign '±'": 1631,
        "Superscript power 's²'": 1637,
        "Newton's dot notation 'ṡ'": 1665,
        "Leibniz's derivative notation 'd/ds'": 1684,
        "At symbol '@' on typewriters": 1884
    }

    print("Symbol introduction years:")
    for symbol, year in symbol_years.items():
        print(f"- {symbol}: {year}")

    # The first possible year is the latest of these dates.
    first_possible_year = max(symbol_years.values())
    print(f"\nThe latest symbol introduction determines the earliest possible year: {first_possible_year}")

    # Round the year to the nearest 10.
    # We use the formula int(year / 10 + 0.5) * 10 for standard rounding.
    rounded_year = int(first_possible_year / 10 + 0.5) * 10
    
    print(f"Rounding {first_possible_year} to the nearest 10 years gives: {rounded_year}")
    
    # The final answer format requires printing the numbers in the final equation.
    # The final "equation" is the rounding calculation.
    print(f"\nFinal calculation: int({first_possible_year} / 10 + 0.5) * 10 = {rounded_year}")

solve_equation_year()