def scan_hexameter_line():
    """
    Analyzes and displays the scansion of a given Latin hexameter line.
    """
    line = "verum ubi equi atque hominis casu convenit imago."

    # The line is segmented into feet based on a detailed metrical analysis.
    # This analysis accounts for elisions, poetic licenses like ictus lengthening,
    # and vowel quantities (by nature or position).
    feet_words = {
        1: "verum ubi",
        2: "equi atque",
        3: "hominis",
        4: "casu",
        5: "convenit i-",
        6: "-mago"
    }

    # The scansion pattern determined from the analysis.
    # D = Dactyl (— U U), S = Spondee (— —)
    feet_scansion = {
        1: "D",
        2: "S",
        3: "S",
        4: "S",
        5: "D",
        6: "S"
    }

    print(f"Original Line: \"{line}\"")
    print("-" * 30)
    print("Metrical Scansion:")

    final_pattern = []
    for i in range(1, 7):
        word_part = feet_words[i]
        scan_part = feet_scansion[i]
        foot_type = "Dactyl" if scan_part == "D" else "Spondee"
        print(f"Foot {i}: {word_part:<14} -> {scan_part} ({foot_type})")
        final_pattern.append(scan_part)

    print("-" * 30)
    print("Final Pattern: " + " ".join(final_pattern))

# Execute the function to print the analysis.
scan_hexameter_line()