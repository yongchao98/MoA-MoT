def find_piazza_feature():
    """
    This function provides information about a historical feature
    in Rome's Piazza della Rotonda.
    """
    feature = "A circular tram line"
    removal_year = 1950
    description = (
        "Until the year {}, the Piazza della Rotonda was encircled by a distinctive {}."
        .format(removal_year, feature.lower())
    )
    details = (
        "These tram tracks looped around the central fountain (Fontana del Pantheon) "
        "and were a prominent part of the square's landscape before being "
        "dismantled to reduce traffic and preserve the historic area."
    )

    print("The unique architectural feature removed from Piazza della Rotonda was:")
    print(f"- Feature: {feature}")
    print(f"- Year of Removal: Approximately {removal_year}")
    print(f"- Description: {description}")
    print(f"- Details: {details}")

find_piazza_feature()