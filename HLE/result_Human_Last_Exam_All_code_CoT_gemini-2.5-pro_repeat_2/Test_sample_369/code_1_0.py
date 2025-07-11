def scan_hexameter_line():
    """
    This function prints the scansion of a specific Latin hexameter line.
    The scansion has been pre-determined due to the line's irregularity and
    the complex rules required to scan it.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The scansion is determined through classical scansion rules,
    # including elisions and some irregular treatments (e.g., for 'equi').
    # D = Dactyl (Long, Short, Short)
    # S = Spondee (Long, Long)
    feet = ["D", "S", "D", "S", "D", "S"]
    
    print(f"Latin line: \"{line}\"")
    print("Scansion (Dactyls as 'D', Spondees as 'S'):")
    
    # Joining the feet with spaces for clarity
    print(" ".join(feet))

scan_hexameter_line()
<<<D S D S D S>>>