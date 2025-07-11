import math

def solve_year_puzzle():
    """
    Determines the earliest possible year an equation could have been written
    based on the notations and typesetting technology used.
    """
    # Step 1: Define the introduction years for key technologies/notations.
    # While other symbols exist (± from 1631, exponents from 1637), they are older
    # and not the limiting factor.
    
    technologies = {
        "Leibniz's derivative notation (d/ds)": 1684,
        "Newton's dot notation (r-dot, s-dot) in print": 1690,
        "TeX typesetting system (Computer Modern font)": 1978
    }

    print("Step-by-step reasoning:")
    print("1. Identifying the key notations and technologies and their introduction years:")
    for tech, year in technologies.items():
        print(f"   - {tech}: {year}")

    # Step 2: The earliest possible year is determined by the most recent invention.
    earliest_possible_year = max(technologies.values())
    
    limiting_tech = max(technologies, key=technologies.get)
    print(f"\n2. The limiting factor is the most recent technology: '{limiting_tech}'.")
    print(f"   This sets the earliest possible year to {earliest_possible_year}.")

    # Step 3: Round the year to the nearest 10 years as requested.
    # We use a custom rounding function for clarity. A positive number rounds up at 5.
    final_year = round(earliest_possible_year / 10) * 10
    
    print("\n3. Rounding this year to the nearest 10 years.")
    
    # Step 4: Display the final calculation with all numbers.
    print("\nFinal Calculation:")
    print(f"The equation uses the font and typesetting style of TeX, invented in {technologies['TeX typesetting system (Computer Modern font)']}.")
    print(f"This is the latest of all technologies present.")
    print(f"Rounding {earliest_possible_year} to the nearest 10 gives:")
    # Using math.ceil for rounding .5 up, which is a common standard.
    # int(earliest_possible_year / 10 + 0.5) * 10 would also work for positive numbers.
    print(f"{final_year} = round({earliest_possible_year} to the nearest 10)")

solve_year_puzzle()
<<<1980>>>