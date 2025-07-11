def scan_hexameter():
    """
    Provides the scansion for the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    
    The scansion is represented by 'D' for a dactyl (— U U) and 'S' for a spondee (— —).
    """
    
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # The established scansion for this complex line is:
    # Dactyl - Spondee - Spondee - Spondee - Dactyl - Spondee
    scansion_feet = ["D", "S", "S", "S", "D", "S"]
    
    print(f"Line: \"{line}\"")
    print("Scansion (D = Dactyl, S = Spondee):")
    
    # The prompt asks to "output each number in the final equation".
    # I will interpret this as outputting each foot's mark in a sequence.
    foot1 = scansion_feet[0]
    foot2 = scansion_feet[1]
    foot3 = scansion_feet[2]
    foot4 = scansion_feet[3]
    foot5 = scansion_feet[4]
    foot6 = scansion_feet[5]
    
    print(f"{foot1} + {foot2} + {foot3} + {foot4} + {foot5} + {foot6}")

# Execute the function to print the result.
scan_hexameter()