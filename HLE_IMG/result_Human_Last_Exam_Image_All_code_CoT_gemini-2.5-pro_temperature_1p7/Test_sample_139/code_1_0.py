import math

def solve_year_puzzle():
    """
    Calculates the earliest possible year for the equation based on its symbols.
    """
    # A dictionary of mathematical notations and their first year of publication.
    notations = {
        'Equals sign "="': 1557,
        'Plus-minus sign "±"': 1631,
        'Exponent notation "x²"': 1637,
        'Leibniz derivative notation "d/dx"': 1684,
        'Newton dot notation "ẋ"': 1736
    }

    # The first possible year is limited by the most recently introduced notation.
    first_possible_year = 0
    limiting_notation = ""
    for notation, year in notations.items():
        if year > first_possible_year:
            first_possible_year = year
            limiting_notation = notation

    # The problem requires rounding to the nearest 10 years.
    # The formula round(number, -1) rounds a number to the nearest 10 in Python 3.
    rounded_year = int(round(first_possible_year, -1))

    print("Analyzing the symbols in the equation to find their introduction year:")
    for notation, year in notations.items():
      print(f"- The {notation} was first published around {year}.")
    
    print(f"\nThe equation could only be written after all its symbols were known.")
    print(f"The limiting factor is the {limiting_notation}, introduced in {first_possible_year}.")
    print(f"Rounding {first_possible_year} to the nearest 10 years gives {rounded_year}.")

solve_year_puzzle()
