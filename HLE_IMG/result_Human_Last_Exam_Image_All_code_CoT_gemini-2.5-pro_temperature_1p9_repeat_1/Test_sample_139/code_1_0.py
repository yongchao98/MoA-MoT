import math

def find_earliest_publication_year():
    """
    Calculates the earliest possible year an equation could be written based on its symbols.
    """
    # A dictionary mapping mathematical symbols/notations to their year of introduction or widespread availability.
    symbol_introduction_years = {
        'Equals sign (=)': 1557,
        'Plus-minus sign (±)': 1631,
        'Superscript exponent notation (x²)' : 1637,
        'Leibniz notation for derivatives (d/dx)': 1684,
        'Newton notation for derivatives (ẋ)': 1693,
        'At sign (@) on a commercial typewriter': 1882
    }

    print("Historical introduction years for symbols in the equation:")
    for symbol, year in symbol_introduction_years.items():
        print(f"- {symbol}: {year}")

    # The equation could not have been written before the last symbol was introduced.
    # Find the maximum year from the dictionary values.
    latest_year = max(symbol_introduction_years.values())

    print(f"\nThe latest introduction year among all symbols is {latest_year}.")
    print("This is the first possible year the equation could have been written.")

    # Round the year to the nearest 10 years.
    # The formula is: round(number / 10) * 10
    rounding_divisor = 10
    rounding_multiplier = 10
    rounded_year = round(latest_year / rounding_divisor) * rounding_multiplier

    print(f"\nRounding {latest_year} to the nearest 10 years.")
    # The final instruction asks to output each number in the final equation.
    # The equation is: rounded_year = round(latest_year / 10) * 10
    print(f"The calculation is: round({latest_year} / {rounding_divisor}) * {rounding_multiplier} = {rounded_year}")

if __name__ == "__main__":
    find_earliest_publication_year()
    # The final rounded year is 1880.
    # Providing the final answer in the requested format.
    print("\n<<<1880>>>")
