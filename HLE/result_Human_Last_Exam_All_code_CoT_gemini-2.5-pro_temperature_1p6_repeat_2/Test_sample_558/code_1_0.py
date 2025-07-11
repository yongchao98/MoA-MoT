def solve_opera_riddle():
    """
    This script solves a multi-step opera history riddle.
    It identifies the opera, the specific production, and the singer in a key role.
    """

    # Step 1: Identify the opera based on the Marietta Alboni clue.
    alboni_year = 1843
    opera_revived_by_alboni = "Semiramide"
    print(f"In {alboni_year}, Marietta Alboni famously sang in a revival of Rossini's opera '{opera_revived_by_alboni}' at La Scala.")

    # Step 2: Establish the timeline based on the Enrico Caruso clue.
    caruso_last_met_year = 1920
    time_gap_years = 70
    revival_target_year = caruso_last_met_year + time_gap_years
    print(f"Enrico Caruso's last Met performance was in {caruso_last_met_year}.")
    print(f"The clue specifies a production more than {time_gap_years} years later, meaning after {revival_target_year}.")

    # Step 3: Identify the specific NYC production.
    nyc_revival_year = 1990
    print(f"The Metropolitan Opera revived '{opera_revived_by_alboni}' in {nyc_revival_year}, which fits the timeframe.")

    # Step 4: Identify the singer of the bass role in that production.
    # The principal bass role in 'Semiramide' is Assur.
    bass_singer = "Samuel Ramey"
    print(f"The principal bass role in the {nyc_revival_year} production was sung by:")
    print(f"\nAnswer: {bass_singer}")

solve_opera_riddle()
<<<Samuel Ramey>>>