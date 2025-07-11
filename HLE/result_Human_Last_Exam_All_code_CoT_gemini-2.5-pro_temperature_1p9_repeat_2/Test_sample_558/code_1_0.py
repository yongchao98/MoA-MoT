def solve_opera_riddle():
    """
    This function solves the riddle by simulating the research steps
    to find the name of the bass singer.
    """
    # Step 1: Identify the opera based on the historical clues provided.
    # Marietta Alboni, La Scala, 1843 -> Rossini's "Semiramide".
    opera = "Semiramide"
    
    # Step 2: Identify the specific NYC production.
    # The Met's 1990 revival fits the timeframe "more than 70 years after Caruso's Met career".
    nyc_production = {
        "venue": "Metropolitan Opera",
        "year": 1990,
        "opera": opera
    }

    # Step 3: Identify the singer of the principal bass role (Assur) in that production.
    # A search for the cast list of the 1990 Met production provides the answer.
    production_cast = {
        "bass_role": "Assur",
        "bass_singer": "Samuel Ramey"
    }

    print(f"The opera revived for Marietta Alboni at La Scala in 1843 was {opera}.")
    print(f"The New York City production referred to is the {nyc_production['year']} revival at the {nyc_production['venue']}.")
    print(f"The principal bass role in this opera is {production_cast['bass_role']}.")
    print(f"In the {nyc_production['year']} production, this role was sung by: {production_cast['bass_singer']}")

solve_opera_riddle()