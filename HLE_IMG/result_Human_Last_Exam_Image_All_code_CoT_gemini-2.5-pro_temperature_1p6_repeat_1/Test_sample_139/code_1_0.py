import math

def main():
    """
    This script determines the earliest possible year an equation could have been written
    based on the history of its symbols, and rounds the year to the nearest 10.
    """

    # Dates of introduction/adoption for the symbols in the equation.
    symbol_dates = {
        "Equals sign (=)": 1557,
        "Plus-minus sign (±)": 1631,
        "Calculus notation (d/dx, dot notation)": 1675,
        "At symbol (@) on commercial typewriters": 1885
    }

    print("Analyzing the history of the symbols:")
    for symbol, year in symbol_dates.items():
        print(f"- {symbol}: ~{year}")

    # The first possible year is determined by the most recently introduced symbol.
    first_possible_year = max(symbol_dates.values())
    print(f"\nThe latest symbol adoption dictates the earliest possible year for the equation.")
    print(f"First possible year: {first_possible_year}")

    # Round the year to the nearest 10.
    # We add 5 and then use integer division by 10, then multiply by 10.
    # This correctly rounds a year like 1885 up to 1890.
    rounding_base = 10
    rounding_add = rounding_base // 2
    rounded_year = (first_possible_year + rounding_add) // rounding_base * rounding_base

    print(f"Rounding {first_possible_year} to the nearest {rounding_base} years...")
    print(f"Calculation: ({first_possible_year} + {rounding_add}) // {rounding_base} * {rounding_base} = {rounded_year}")
    print(f"\nThe final answer is: {rounded_year}")

if __name__ == "__main__":
    main()
