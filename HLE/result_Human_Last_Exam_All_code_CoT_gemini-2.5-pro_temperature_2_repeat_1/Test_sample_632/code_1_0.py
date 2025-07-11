def solve_trivia():
    """
    This function solves the trivia question about the Piazza della Rotonda.
    """
    feature = "The tramway lines"
    removal_year = 1950
    explanation = (
        f"Until around {removal_year}, {feature.lower()} ran directly through the Piazza della Rotonda, "
        "passing in front of the Pantheon. These tracks and their overhead power lines were a "
        "defining feature of the square's modern landscape for decades. They were removed as "
        "part of a city-wide effort to update transportation and reduce the visual impact of "
        "modern infrastructure on ancient monuments."
    )
    
    print("The unique feature removed from the Piazza della Rotonda around 1950 was:")
    print(explanation)

solve_trivia()