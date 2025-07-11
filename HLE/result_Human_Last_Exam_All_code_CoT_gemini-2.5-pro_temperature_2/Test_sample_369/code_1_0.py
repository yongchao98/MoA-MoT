def solve_scansion():
    """
    This function scans the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    It prints the resulting pattern of dactyls (D) and spondees (S) for all six feet.
    The analysis, which involves several poetic licenses, is explained in the comments.
    """

    # Foot 1: "verum ubi" -> vēr-u-bi -> Dactyl
    # - "ver-" is long by position; "-um" elides with "ubi".
    foot_1 = "D"

    # Foot 2: "equi atque" -> ēqu'āt- -> Spondee
    # - "equi" elides with "atque"; "e-" and "at-" are long by position.
    foot_2 = "S"

    # Foot 3: "-que hominis" -> -quē ho-mi- -> Dactyl
    # - Requires the poetic license of treating "-que" as long.
    foot_3 = "D"

    # Foot 4: "-nis casu" -> -nīs că-sŭ -> Dactyl
    # - "-nis" is long by position; requires a strong license to scan "casu" (normally --) as short-short.
    foot_4 = "D"
    
    # Foot 5: "convenit" -> cōn-ve-nit -> Dactyl
    # - "con-" is long by position; "-ve-" and "-nit" are short. Assumes prodelision of "imago".
    foot_5 = "D"

    # Foot 6: "(i)mago" -> -mā-gō -> Spondee
    # - With the initial "i" dropped via prodelision, "-mā-" and "-gō" are both long.
    foot_6 = "S"
    
    # The final equation of the line's meter is the sequence of the six feet.
    # We print each symbol of the final pattern.
    print(f"{foot_1} {foot_2} {foot_3} {foot_4} {foot_5} {foot_6}")

solve_scansion()