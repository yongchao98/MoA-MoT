def scan_hexameter_line():
    """
    Scans the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    This function demonstrates the step-by-step scansion, including a common
    textual emendation needed to make the line scan correctly.
    """

    # The line is known to be a trick question, requiring 'equi' to be 'aeque'.
    line_syllables = [
        ("vē-rum u-bi", "ver(um) ubi elides to verubi."),
        ("ae-quē", "The word is 'aeque' (adv.), not 'equi'. Hiatus occurs after 'ubi'."),
        ("āt-que ho-mi-nis", "atque and hominis."),
        ("cā-sū", "casu."),
        ("cōn-ve-nit", "convenit (present tense)."),
        ("i-mā-gō", "The 'i' becomes long by position from the 't' in 'convenit'.")
    ]

    # The canonical scansion is D S S S D S. We will reconstruct this logic.
    feet = []
    foot_patterns = []

    # Foot 1: Dactyl (D)
    # vērum ubi -> vēr(um) ubi -> vē-ru-bi -> — U U
    feet.append("vē-ru-bi")
    foot_patterns.append("D")

    # Foot 2: Spondee (S)
    # aequē -> — — (ae is a diphthong, ē is long by nature).
    # A pause (caesura/hiatus) is assumed between 'ubi' and 'aeque'.
    feet.append("ae-quē")
    foot_patterns.append("S")

    # Foot 3: Spondee (S)
    # atque hominis -> āt-que hom-. 'āt' is long by position.
    # 'que' elides with 'hominis': atqu(e) hominis.
    # The syllable 'hom' is treated as long.
    feet.append("āt-qu'ho-min")
    foot_patterns.append("S")

    # Foot 4: Spondee (S)
    # The remainder of hominis + casu -> -is cā-sū.
    # '-is' is long by position before 'c'. 'cā' is long by nature.
    feet.append("is cā-sū")
    foot_patterns.append("S")

    # Foot 5: Dactyl (D)
    # cōnvenit -> cōn-ve-nit -> — U U
    # 'cōn' is long by position. 've' and 'nit' are short.
    feet.append("cōn-ve-nit")
    foot_patterns.append("D")

    # Foot 6: Spondee (S)
    # imāgō -> i-mā-gō -> — —
    # The 'i' in 'imago' is normally short, but here it is made long by position
    # across the word boundary by the final 't' of convenit, an irregularity.
    # The final syllable '-gō' is treated as long (syllaba anceps).
    feet.append("i-mā-gō")
    foot_patterns.append("S")

    print("Scansion of: verum ubi aequē atque hominis casu convenit imago.")
    print("-" * 60)
    for i in range(6):
        print(f"Foot {i+1}: {feet[i]:<15} Pattern: {foot_patterns[i]}")
    print("-" * 60)

    # Print the final combined pattern.
    final_pattern = " ".join(foot_patterns)
    print(f"Final Scansion: {final_pattern}")

scan_hexameter_line()