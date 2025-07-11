import math

def solve_year_puzzle():
    """
    Determines the earliest possible year an equation could have been written
    based on the introduction dates of its symbols.
    """
    # Historical dates for the introduction of mathematical symbols and notations.
    # While most symbols are from the 16th-17th centuries, the use of '@' as a
    # variable in a font that looks like Computer Modern is the key. This points
    # to the TeX typesetting system, first released in 1978.
    symbol_introduction_years = {
        "Equals sign '='": 1557,
        "Plus-minus sign '±'": 1631,
        "Exponent notation 's²'": 1637,
        "Newton's dot notation 'ṡ'": 1665,
        "Leibniz's differential 'd/ds'": 1675,
        "At-sign '@' in TeX typesetting": 1978,
    }

    # The first possible year the equation could have been written is determined
    # by the latest of these dates.
    first_possible_year = max(symbol_introduction_years.values())

    # The problem asks to round the year to the nearest 10 years.
    # In Python, rounding to the nearest 10 can be done by dividing by 10,
    # rounding to the nearest integer, and then multiplying by 10.
    base = 10
    rounded_year = int(base * round(float(first_possible_year) / base))

    print(f"The latest symbol introduction date is: {first_possible_year}")
    print(f"The final calculation is rounding the year {first_possible_year} to the nearest {base} years.")
    print(f"Result: {rounded_year}")

solve_year_puzzle()