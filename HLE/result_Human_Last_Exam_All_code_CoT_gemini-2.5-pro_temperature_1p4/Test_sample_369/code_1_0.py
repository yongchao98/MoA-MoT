def scan_hexameter():
    """
    Prints the scansion of the given Latin hexameter line.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    scansion = ["D", "S", "D", "S", "D", "S"]
    
    print(f"Latin line: \"{line}\"")
    print("Scansion (D=Dactyl, S=Spondee):")
    
    # Joining with ' - ' to clearly separate the feet
    result = " - ".join(scansion)
    print(result)

scan_hexameter()