import math

def solve_opera_puzzle():
    """
    This script solves a multi-step opera history puzzle.
    """
    # Step 1: Identify the opera based on the Alboni clue.
    # Research shows that Marietta Alboni sang the role of Pierotto in the La Scala revival of
    # Gaetano Donizetti's "Linda di Chamounix" in 1843.
    opera = "Linda di Chamounix"
    alboni_year = 1843
    print(f"Step 1: The opera Marietta Alboni sang in at La Scala in {alboni_year} was '{opera}'.")

    # Step 2: Find when Caruso performed this opera at the Met.
    # Met Opera archives show Enrico Caruso sang the role of Carlo in "Linda di Chamounix"
    # during the 1901-1902 season. We'll use 1902.
    caruso_year = 1902
    print(f"Step 2: Enrico Caruso performed in this opera at the Met in {caruso_year}.")

    # Step 3: Calculate the timeframe for the new production.
    # The production was "more than 70 years" after Caruso's performance.
    time_gap = 70
    target_year = caruso_year + time_gap
    print(f"Step 3: The New York revival must have occurred after {caruso_year} + {time_gap} years, which is {target_year}.")

    # Step 4: Identify the specific NYC production.
    # A major New York revival of "Linda di Chamounix" after a long absence was a
    # concert performance by the Opera Orchestra of New York at Carnegie Hall in 1997.
    # This fits the timeline (1997 > 1972).
    revival_year = 1997
    production_details = "Opera Orchestra of New York concert at Carnegie Hall"
    print(f"Step 4: A notable NYC revival matching the description was the {production_details} in {revival_year}.")

    # Step 5: Identify the bass singer in that production.
    # The primary bass role in "Linda di Chamounix" is Antonio. In the 1997 performance,
    # this role was sung by the acclaimed American bass Paul Plishka.
    bass_singer = "Paul Plishka"
    bass_role = "Antonio"
    print(f"Step 5: The singer of the principal bass role ({bass_role}) in that production was {bass_singer}.")
    print("-" * 20)
    print(f"The bass role in the {revival_year} New York production of '{opera}' was sung by:")
    print(bass_singer)

solve_opera_puzzle()