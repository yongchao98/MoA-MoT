def solve_riddle():
    """
    This script solves the riddle by breaking it down into its components.
    """

    # Step 1: Identify Carl Sagan's "X of Y" phrase.
    # The riddle describes Sagan's view on Earth's images from space.
    # The most fitting phrase is "Shore of the Cosmic Ocean".
    X = "Shore"
    Y = "Ocean"

    # Step 2: Identify the crater and its name's etymology.
    # The crater in the image is the "Sagan" crater on Mars.
    # The riddle connects the name's etymology to the phrase.
    # A plausible Old English origin for the name "Sagan" links it to "sæ", meaning "sea" or "ocean".
    crater_name = "Sagan"
    etymological_root_meaning = "sea/ocean"

    # Step 3: Conclude the answer.
    # The phrase is "Shore of the Ocean". Therefore, Y is "Ocean".
    print(f"Carl Sagan's famous phrase, fitting the 'X of Y' format, is '{X} of the cosmic {Y}'.")
    print(f"The crater in the image is named '{crater_name}'.")
    print(f"The name '{crater_name}' has a plausible etymological link to the Old English word for '{etymological_root_meaning}'.")
    print(f"This connects the man to his metaphor.")
    print("\nTherefore, Y is:")
    print(Y)

solve_riddle()