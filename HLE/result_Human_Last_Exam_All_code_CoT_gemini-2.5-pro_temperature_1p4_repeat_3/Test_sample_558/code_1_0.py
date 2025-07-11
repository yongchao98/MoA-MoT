import sys

def solve_opera_riddle():
    """
    This function solves the opera riddle by following a logical deduction process.
    """
    # Step 1: Identify the opera based on the Alboni clue.
    opera = "Semiramide"
    year_alboni_la_scala = 1843
    print(f"Step 1: The opera revived at La Scala in {year_alboni_la_scala} with Marietta Alboni was Rossini's '{opera}'.")

    # Step 2: Establish a timeline using the Caruso clue.
    # Enrico Caruso's last full season at the Met was 1920.
    caruso_era_end = 1920
    time_gap = 70
    production_year = caruso_era_end + time_gap
    print(f"Step 2: The production was staged more than {time_gap} years after Caruso's era.")
    print("Using his last full season at the Met as a reference point.")
    
    # Step 3: Print the calculation as requested.
    print("\nStep 3: Calculating the target year for the New York revival:")
    print(f"Calculation: {caruso_era_end} + {time_gap} = {production_year}")
    print(f"This points to a significant production of '{opera}' around or after {production_year}.")

    # Step 4: Identify the specific production and the bass singer.
    # A famous revival of "Semiramide" premiered at the Metropolitan Opera in 1990.
    # The principal bass role in this opera is Assur.
    nyc_production_year = 1990
    bass_role = "Assur"
    bass_singer = "Samuel Ramey"
    
    print(f"\nStep 4: The Metropolitan Opera staged a landmark revival of '{opera}' in {nyc_production_year}.")
    print(f"The singer who performed the bass role of '{bass_role}' in this production was {bass_singer}.")
    
    # Final Answer
    print("\n---")
    print(f"Final Answer: The bass role in the {nyc_production_year} New York City production of '{opera}' was sung by:")
    print(bass_singer)

# Execute the function
if __name__ == "__main__":
    solve_opera_riddle()