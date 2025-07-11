def find_earliest_publication_year():
    """
    Calculates the earliest possible year an equation could have been written based on the history of its symbols.
    """
    # Step 1 & 2: A dictionary mapping symbols/notations to their year of introduction.
    # The key is the availability of the symbol for printing in a published paper.
    symbol_years = {
        "Equals sign (=)": 1557,
        "Plus-Minus sign (±)": 1631,
        "Newton's notation (dot)": 1665,
        "Leibniz's notation (d/ds)": 1684,
        "'@' symbol on a commercial typewriter": 1874
    }

    print("To determine the earliest possible year, we analyze the introduction year of each symbol in the equation:")
    print("-" * 40)
    for symbol, year in symbol_years.items():
        print(f"{symbol}: {year}")
    print("-" * 40)

    # Step 3: Find the latest year among the symbols.
    # The paper could not have been published before the last symbol was available.
    latest_year = 0
    all_years = []
    for year in symbol_years.values():
        all_years.append(year)
        if year > latest_year:
            latest_year = year
            
    # Showing the numbers used in the max calculation
    print(f"The first possible year is the latest of these dates.")
    print(f"Calculation: max({', '.join(map(str, all_years))})")
    print(f"The latest year is: {latest_year}")

    # Step 4: Round the year to the nearest 10 years.
    # round(n / 10.0) * 10 is a standard way to do this.
    rounded_year = int(round(latest_year / 10.0) * 10)
    
    print("\nThe problem requires rounding this year to the nearest 10 years.")
    print(f"Rounding calculation: round({latest_year} / 10) * 10 = round({latest_year / 10.0}) * 10 = {int(round(latest_year / 10.0))} * 10 = {rounded_year}")
    
    return rounded_year

# Execute the function and print the final answer.
final_answer = find_earliest_publication_year()
print(f"\nThe first possible year this equation could have been written, rounded to the nearest 10 years, is {final_answer}.")
print(f"\n<<<1870>>>")
