def scan_lucretius_line():
    """
    Prints the scansion of the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    This specific line requires several poetic licenses to be scanned properly.
    The final scansion is represented by 'D' for dactyls and 'S' for spondees.
    """
    line = {
        "verum ubi": "D",
        "equi atque": "S",
        "hominis": "D",
        "casu": "S",
        "convenit": "D",
        "imago": "S"
    }

    scansion_parts = []
    print("Scansion of the line:")
    print("verum ubi equi atque hominis casu convenit imago.")
    print("-" * 30)
    for part, foot_type in line.items():
        # The feet do not align perfectly with these word groups,
        # but this representation is clear.
        # For example, the final spondee is formed by syllables from both convenit and imago.
        print(f"'{part}' scans as: {foot_type}")
        scansion_parts.append(foot_type)

    final_pattern = " ".join(scansion_parts)
    print("-" * 30)
    print("The final pattern for the six feet is:")
    print(final_pattern)

scan_lucretius_line()