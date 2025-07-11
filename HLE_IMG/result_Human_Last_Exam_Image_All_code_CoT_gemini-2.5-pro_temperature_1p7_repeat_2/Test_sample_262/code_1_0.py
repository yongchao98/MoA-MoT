def find_shorter_lifetime_complexes():
    """
    Identifies iridium complexes expected to have shorter lifetimes based on
    the position of fluorine substituents.
    """

    # Data representing the fluorine atom positions on the phenyl ring for each complex.
    # The C-Ir bond is at position C2'. The critical 'ortho' position is C3'.
    # A list of numbers represents the positions of the fluorine atoms.
    complex_structures = {
        '1': [4, 5],  # No ortho-fluorine
        '2': [4],     # No ortho-fluorine
        '3': [3],     # Contains an ortho-fluorine
        '4': [3, 4, 5] # Contains an ortho-fluorine
    }

    destabilizing_position = 3
    shorter_lifetime_list = []

    print("Analyzing complexes based on the rule: Fluorine at the ortho-position (C3') shortens lifetime.")
    print("-" * 30)

    # Sort keys to ensure consistent output order
    for complex_id in sorted(complex_structures.keys()):
        positions = complex_structures[complex_id]
        if destabilizing_position in positions:
            shorter_lifetime_list.append(int(complex_id))
            print(f"Complex {complex_id}: Found a destabilizing fluorine at position {destabilizing_position}. -> Shorter lifetime expected.")
        else:
            print(f"Complex {complex_id}: No fluorine at position {destabilizing_position}. -> Longer lifetime expected.")

    print("-" * 30)
    print("Conclusion: The complexes expected to show shorter lifetimes are:")
    # Using print to show the final equation/list as requested.
    # We output each number in the final selection.
    print(f"Final List = [{', '.join(map(str, shorter_lifetime_list))}]")


find_shorter_lifetime_complexes()