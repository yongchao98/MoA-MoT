def scan_hexameter_line():
    """
    Prints the scansion of a given Latin hexameter line.
    The scansion represents dactyls as 'D' and spondees as 'S'.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # The scansion for this line is Dactyl-Dactyl-Spondee-Spondee-Spondee-Spondee
    scansion_feet = ["D", "D", "S", "S", "S", "S"]
    
    print(f"Latin line: {line}")
    
    # Print the final scansion, showing each foot's type
    final_scansion = " - ".join(scansion_feet)
    print(f"Scansion: {final_scansion}")

scan_hexameter_line()