def find_borges_names():
    """
    This function identifies and prints the three names mentioned by Jorge Luis Borges
    in his essay "Our Poor Individualism" to illustrate the illusions of patriotism.
    """
    # The three individuals cited by Borges in the essay.
    greek_figure = "Plutarch"
    english_figure = "Gibbon"
    germanist_figure = "Thorndike"

    # Combine the names into a single string, separated by commas.
    names = [greek_figure, english_figure, germanist_figure]
    result = ", ".join(names)

    # Print the final result.
    print(result)

find_borges_names()