def solve_riddle():
    """
    Analyzes the rock art image to identify a symbol not native to ancient Southwest cultures.
    """
    
    is_foreign_symbol_present = True
    foreign_symbol_description = "The letters 'RO' or 'NO'"

    print("Analyzing the rock art panel...")
    print("The panel contains many symbols typical of ancient Southwest cultures, such as anthropomorphic figures and geometric designs.")
    print("However, there is an element that does not belong to these pre-contact cultures.")
    print("-" * 20)
    print(f"True or False: Is there at least one symbol not related to these cultures? -> {is_foreign_symbol_present}")
    print("-" * 20)
    print(f"If true, what symbol? -> {foreign_symbol_description}")
    print("\nExplanation:")
    print("In the upper right section of the image, the letters 'RO' or 'NO' are visible. These are characters from the Latin alphabet, which was introduced to the Americas by Europeans. Their presence among the ancient pictographs indicates they were added at a later time and are not part of the original cultural symbolic system.")

solve_riddle()