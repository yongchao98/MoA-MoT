def solve_opera_riddle():
    """
    This script solves a multi-step opera trivia question.
    """
    # Step 1: Identify the opera based on the Marietta Alboni clue.
    opera = "Rossini's Semiramide"
    alboni_performance_year = 1843

    # Step 2: Identify the specific NYC production based on the timeline.
    # The Met last performed Semiramide in 1895.
    last_met_performance_year = 1895
    # The next major NYC revival was at the Met in 1990.
    revival_year = 1990
    
    # The time gap fits the clue "more than 70 years".
    time_gap = revival_year - last_met_performance_year
    
    # Step 3: Identify the bass singer in that 1990 production.
    # The principal bass role is Assur.
    bass_singer = "Samuel Ramey"

    # Print the final answer
    print(f"The opera is {opera}.")
    print(f"The revival production in New York City was at the Metropolitan Opera in {revival_year}.")
    print(f"The time between the previous Met production ({last_met_performance_year}) and this revival was {time_gap} years.")
    print(f"The bass role in this production was sung by: {bass_singer}")

solve_opera_riddle()