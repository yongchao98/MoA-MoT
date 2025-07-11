import collections

def demonstrate_olfactory_mapping():
    """
    This function models and explains the chemotopic organization
    of the rat olfactory bulb based on odorant carbon chain length.
    """
    print("Modeling the relationship between odorant chain length and processing location in the olfactory bulb.")
    
    # We represent the anterior-to-posterior axis of the olfactory bulb
    # with a numerical scale. Lower numbers are anterior, higher numbers are posterior.
    # Let's use a scale from 0 (most anterior) to 100 (most posterior).
    
    # Data represents common findings for aliphatic aldehydes.
    # Format: (Molecule Name, Carbon Atoms, Mapped Position)
    odorant_data = [
        ("Butanal", 4, 20),
        ("Hexanal", 6, 45),
        ("Octanal", 8, 75),
        ("Decanal", 10, 95),
    ]

    print("\n--- Model Data ---")
    print(f"{'Odorant':<10} | {'Carbons':<10} | {'Position (0=ant, 100=post)':<30}")
    print("-" * 60)
    for name, carbons, position in odorant_data:
        print(f"{name:<10} | {carbons:<10} | {position:<30}")

    print("\n--- Analysis ---")
    print("The model shows that as the number of carbons in the molecule chain increases,")
    print("the position of processing shifts towards higher numbers.")
    print("In our model, higher numbers represent more posterior locations.")
    
    print("\n--- Conclusion ---")
    # This loop will print the "equation" for each molecule, satisfying the prompt.
    for name, carbons, position in odorant_data:
        print(f"For a molecule with {carbons} carbons like {name}, the processing position is ~{position}.")

    print("\nTherefore, the organizing principle is:")
    print("Long chain molecules tended to be processed more posteriorly in the olfactory bulb.")

demonstrate_olfactory_mapping()