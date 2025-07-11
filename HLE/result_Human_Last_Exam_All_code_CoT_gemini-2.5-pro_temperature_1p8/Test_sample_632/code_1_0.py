def solve_historical_puzzle():
    """
    This function provides the answer to the historical question about
    the Piazza della Rotonda.
    """

    # The year the feature was removed.
    year_removed = 1950

    # The architectural feature in question.
    feature_name = "Tram tracks"

    # A detailed explanation of the feature and its removal.
    explanation = (
        f"Until around the year {year_removed}, a network of tram tracks ran directly through "
        "the Piazza della Rotonda. These tracks, part of Rome's extensive public transit system, "
        "were considered a unique, albeit perhaps noisy, feature of the square's daily life. "
        "They were removed by the city administration as part of a post-war effort to 'modernize' "
        "the historic center and give priority to buses and cars, effectively changing the square's character."
    )

    print("The unique architectural feature removed from Piazza della Rotonda was:")
    print(f"- {feature_name}")
    print("\nMore details:")
    print(explanation)

# Execute the function to print the answer.
solve_historical_puzzle()