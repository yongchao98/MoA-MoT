import calendar

def solve_opera_riddle():
    """
    This script solves a multi-step opera history riddle.
    """
    # Step 1: Identify the opera based on the 1843 La Scala revival with Marietta Alboni.
    opera_title = "L'italiana in Algeri"
    alboni_revival_year = 1843
    print(f"The opera revived at La Scala in {alboni_revival_year} with Marietta Alboni was Rossini's '{opera_title}'.")

    # Step 2: Identify the year of Caruso's last performance of this opera at the Met.
    caruso_last_performance_year = 1919
    print(f"Enrico Caruso's last performance in this opera at the Metropolitan Opera was in {caruso_last_performance_year}.")

    # Step 3: Calculate the year of the new production. It's more than 70 years after Caruso's last performance.
    time_gap_years = 70
    target_year = caruso_last_performance_year + time_gap_years
    print(f"\nThe riddle specifies a production that took place more than {time_gap_years} years later.")
    print(f"The calculation for the earliest possible year is: {caruso_last_performance_year} + {time_gap_years} = {target_year}")
    
    # Step 4: Identify the specific NYC production that fits the timeline.
    # The Met staged a new production for the first time since Caruso's era in 1991.
    nyc_production_year = 1991
    years_since_caruso = nyc_production_year - caruso_last_performance_year
    print(f"The first New York City staging after that long gap was the Metropolitan Opera production in {nyc_production_year}.")
    print(f"This was {years_since_caruso} years after Caruso's last performance, satisfying the 'more than 70 years' condition.")

    # Step 5: Identify the bass singer in that production.
    # The principal bass role in this opera is Mustafà.
    bass_singer = "Samuel Ramey"
    bass_role = "Mustafà"
    print(f"\nThe bass role of {bass_role} in the {nyc_production_year} Met Opera production was sung by Samuel Ramey.")

solve_opera_riddle()
<<<Samuel Ramey>>>