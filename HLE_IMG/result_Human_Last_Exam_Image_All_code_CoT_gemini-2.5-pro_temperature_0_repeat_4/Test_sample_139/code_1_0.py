import math

def find_earliest_publication_year():
    """
    Determines the earliest possible publication year of an equation based on its symbols.
    """
    # A dictionary mapping symbols/notations to their year of first significant appearance or availability.
    symbol_introduction_years = {
        "Equals sign (=)": 1557,
        "Plus-minus sign (±)": 1631,
        "Superscript power notation (x²)": 1637,
        "Leibniz's differential notation (d/dx)": 1684,
        "Newton's dot notation (ẋ)": 1693,
        "At symbol (@) on a commercial typewriter": 1896
    }

    print("Analyzing the symbols in the equation to find the latest one introduced:")
    
    latest_year = 0
    limiting_symbol = ""

    # Print each symbol and its year, and find the maximum year
    for symbol, year in symbol_introduction_years.items():
        print(f"- The {symbol} was introduced around the year {year}.")
        if year > latest_year:
            latest_year = year
            limiting_symbol = symbol

    print(f"\nThe symbol that appeared last is the '{limiting_symbol}'.")
    print(f"Therefore, the earliest the equation could have been written is {latest_year}.")

    # Round the latest year to the nearest 10
    rounded_year = int(round(latest_year / 10.0)) * 10

    print(f"\nRounding {latest_year} to the nearest 10 years results in {rounded_year}.")
    
    # Per instructions, showing the numbers in the final equation
    print("\nFinal Equation:")
    years = list(symbol_introduction_years.values())
    print(f"Year = round(max({', '.join(map(str, years))}) / 10) * 10")
    print(f"Year = round({latest_year} / 10) * 10")
    print(f"Year = {rounded_year}")


find_earliest_publication_year()