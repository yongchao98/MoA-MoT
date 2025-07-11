def solve_riddle():
    """
    This function identifies the thinker based on historical context.
    """
    # Key information from the prompt:
    host = "Léon Bloy"
    artist = "Georges Rouault"
    time_period_start = 1905
    time_period_end = 1909
    description = "A prominent figure of European thought"

    # The historical context points to a specific intellectual circle.
    # During this period, a famous philosopher and his wife were frequent visitors
    # at Léon Bloy's home and were close to Georges Rouault.
    # This philosopher perfectly matches the description.
    thinker = "Jacques Maritain"

    print(f"The writer was Léon Bloy, and the artist was Georges Rouault.")
    print(f"The events described took place between {time_period_start} and {time_period_end}.")
    print(f"The prominent thinker who was part of this circle and wrote the observation was Jacques Maritain.")

solve_riddle()