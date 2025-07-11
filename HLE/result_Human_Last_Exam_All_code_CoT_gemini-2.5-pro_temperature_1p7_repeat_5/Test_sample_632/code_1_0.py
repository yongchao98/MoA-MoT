def solve_trivia():
    """
    This function provides the answer to the historical question about the Piazza della Rotonda.
    """
    feature = "the ring of tram tracks"
    location = "Piazza della Rotonda"
    removed_year = 1950

    answer = (f"The unique feature removed from the {location} in {removed_year} "
              f"was {feature} that circled the central fountain. "
              f"They were removed by the city administration to accommodate the rise of automobile traffic.")

    print(answer)

solve_trivia()