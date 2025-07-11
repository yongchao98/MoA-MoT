def find_removed_feature():
    """
    This function identifies the architectural feature removed from
    Piazza della Rotonda around the year 1950.
    """

    # The year 1950 is a clue. It points to the post-WWII modernization of Rome.
    # While the famous bell towers on the Pantheon were removed, this happened much earlier, in 1883.
    # In the first half of the 20th century, tram lines ran directly through many of Rome's historic squares.
    # After the war, these were largely seen as disruptive and were removed in favor of bus lines.
    feature = "tram tracks"
    year = 1950

    # There is no equation, so we will simply print the components of the answer.
    print(f"Feature: {feature}")
    print(f"Approximate year of removal: {year}")
    print(f"The unique feature that ran through the Piazza della Rotonda until it was removed around {year} was the set of {feature}.")

find_removed_feature()