def solve_opera_mystery():
    """
    This function solves a multi-step trivia question about opera history
    by storing facts in variables and printing the logical conclusion.
    """
    # Step 1: Identify the opera based on the 1843 La Scala revival clue.
    opera_revival_year = 1843
    opera_title = "Semiramide"
    print(f"The opera revived at La Scala in {opera_revival_year} featuring Marietta Alboni was Rossini's '{opera_title}'.")

    # Step 2: Use the timeline to identify the modern NYC production.
    # The last Met production of this opera in the previous era was in 1895.
    last_production_year = 1895
    time_gap = 70
    
    # The new production took place in 1990, which is more than 70 years after 1895.
    nyc_revival_year = 1990

    print(f"A new production was staged in New York in {nyc_revival_year}, which is more than {time_gap} years after the previous Met production in {last_production_year}.")

    # This line prints the numbers involved in the timeline logic, as requested.
    print(f"The historical years in the puzzle are: {opera_revival_year}, {last_production_year}, {time_gap}, and {nyc_revival_year}.")

    # Step 3: Identify the singer of the bass role in the 1990 production.
    # The principal bass role in "Semiramide" is Assur.
    bass_singer = "Samuel Ramey"
    
    print(f"\nThe singer of the principal bass role in the {nyc_revival_year} production was:")
    print(bass_singer)

# Run the function to solve the puzzle.
solve_opera_mystery()