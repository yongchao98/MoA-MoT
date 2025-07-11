def scan_hexameter_line():
    """
    Scans a specific line of Latin hexameter and prints the result.
    The line is: "verum ubi equi atque hominis casu convenit imago."
    This line is known for its metrical complexities, and the scansion
    presented is based on established classical analysis.
    """

    # The original Latin line
    line = "verum ubi equi atque hominis casu convenit imago."

    # The line divided into its six metrical feet
    # Note the splitting of words like 'atque', 'convenit', and 'imago' across feet.
    feet_division = ["verum ubi", "equi at", "que hominis", "casu con", "venit i", "mago"]

    # The scansion pattern for each foot: D for Dactyl, S for Spondee
    scansion_pattern = ["D", "S", "S", "S", "D", "S"]

    # Print the header
    print(f"Scansion of the line: \"{line}\"")
    print("-" * 30)

    # Print each foot with its scansion
    full_scansion_line = []
    for i in range(len(feet_division)):
        foot_text = feet_division[i]
        pattern = scansion_pattern[i]
        # We represent the dactyls and spondees in the final equation
        if pattern == "D":
            foot_representation = "D"
        else: # Spondee
            foot_representation = "S"
        
        full_scansion_line.append(foot_representation)

        print(f"Foot {i+1}: {foot_text:<12} -> {pattern}")
    
    # Print the final scansion on one line
    print("-" * 30)
    print("Final Pattern: " + " | ".join(full_scansion_line))


# Execute the function to print the scansion
scan_hexameter_line()