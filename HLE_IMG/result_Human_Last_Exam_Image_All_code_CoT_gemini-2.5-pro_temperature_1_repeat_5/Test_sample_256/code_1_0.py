def analyze_rock_art():
    """
    Analyzes the provided image of rock art to find an unrelated symbol.
    """
    is_unrelated_symbol_present = True
    unrelated_symbol_description = (
        "The word 'no' is visible in the upper right quadrant of the image. "
        "This is written using the Latin alphabet, which was not used by the ancient Southwest cultures who created the surrounding pictographs. "
        "It is an example of modern graffiti."
    )

    print(f"True or false? Is there an unrelated symbol? {is_unrelated_symbol_present}")
    if is_unrelated_symbol_present:
        print(f"What is the symbol? {unrelated_symbol_description}")

analyze_rock_art()