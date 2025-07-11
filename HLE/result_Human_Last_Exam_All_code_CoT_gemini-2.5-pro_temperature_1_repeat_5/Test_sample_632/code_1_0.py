def solve_pantheon_mystery():
    """
    This function explains the architectural feature removed from the Pantheon
    in the Piazza della Rotonda.
    """
    feature_name = "a pair of bell towers"
    nickname = "'le orecchie d'asino' or 'the donkey's ears'"
    architect = "Bernini (or his school)"
    year_added = "17th century"
    year_removed_actual = 1883
    year_removed_question = 1950

    print(f"The unique architectural feature in the Piazza della Rotonda was {feature_name}.")
    print(f"They were added in the {year_added}, likely by the architect {architect}.")
    print(f"The Romans gave them the derogatory nickname {nickname} because they were considered ugly and out of place on the ancient temple.")
    print("\nIMPORTANT NOTE ON THE DATE:")
    print(f"While the question states the removal year was {year_removed_question}, these widely disliked towers were actually removed much earlier.")
    # The final part of the prompt asks to output each number in the final equation.
    # Since there's no calculation, I will simply print the key year.
    print(f"The correct year of their removal by the city administration was: {year_removed_actual}")

solve_pantheon_mystery()