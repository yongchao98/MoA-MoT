def scan_hexameter_line():
    """
    Scans a difficult line of Latin hexameter and prints the result.

    The line is from Lucan's Pharsalia (1.568):
    "verum ubi equi atque hominis casu convenit imago."

    This line is known for its metrical complexity. The scansion D-S-S-S-D-S is
    the one most commonly cited by scholars, though it requires assumptions
    such as synizesis (treating two vowel sounds as one) and metrical
    lengthening (treating a short syllable as long due to metrical stress).

    The breakdown is as follows:
    - Foot 1: vēr(um) ubi -> Dactyl (— U U)
    - Foot 2: equ(i) atque -> Spondee (— —) (Requires synizesis and elision)
    - Foot 3: hominis -> Spondee (— —) (Requires metrical lengthening of 'ho-')
    - Foot 4: cāsū con- -> Spondee (— —)
    - Foot 5: -vēnit i- -> Dactyl (— U U)
    - Foot 6: -māgō -> Spondee (— —)
    """

    line = "verum ubi equi atque hominis casu convenit imago."
    feet_text = ["verum ubi", "equi atque", "hominis", "casu con-", "-venit i-", "-mago"]
    scansion = ["D", "S", "S", "S", "D", "S"]

    print(f"Scanning line: \"{line}\"\n")

    for i in range(len(scansion)):
        print(f"Foot {i+1} ({feet_text[i]}): {scansion[i]}")

    final_pattern = " ".join(scansion)
    print(f"\nFinal pattern: {final_pattern}")

scan_hexameter_line()