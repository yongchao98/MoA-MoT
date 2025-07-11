def solve_piazza_riddle():
    """
    This function stores and retrieves historical information
    about the Piazza della Rotonda to answer the user's question.
    """
    # Historical data about the feature removed from the Pantheon
    removed_feature_data = {
        'name': "Twin bell towers",
        'nickname': "The donkey's ears (le orecchie d'asino)",
        'architect': "Gian Lorenzo Bernini (common attribution)",
        'reason_for_removal': "Considered out of harmony with the classical facade of the Pantheon.",
        'actual_removal_year': 1883
    }

    feature_name = removed_feature_data['name']
    nickname = removed_feature_data['nickname']
    removal_year = removed_feature_data['actual_removal_year']
    
    print(f"The unique architectural feature in question was a set of {feature_name}.")
    print(f"They were added in the 17th century and were widely known by their derogatory nickname: '{nickname}'.")
    print(f"Although the prompt mentions 1950, historical records show they were actually removed in {removal_year}.")

solve_piazza_riddle()