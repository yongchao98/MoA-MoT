def scan_hexameter_line():
    """
    Prints the scansion for a specific line of Latin hexameter.
    The scansion for this particular line is based on scholarly consensus due to its metrical complexity.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The accepted scansion is Dactyl-Spondee-Spondee-Spondee-Dactyl-Spondee.
    scansion_pattern = ["D", "S", "S", "S", "D", "S"]

    print("Latin Hexameter Line:")
    print(f"'{line}'")
    print("\nScansion (D=Dactyl, S=Spondee):")

    # The prompt requires outputting each symbol in the final result.
    # We will print the full pattern representing all six feet.
    final_output = " - ".join(scansion_pattern)
    print(final_output)

# Execute the function to display the result.
scan_hexameter_line()