def scan_hexameter():
    """
    Scans the provided Latin hexameter line and prints the result.
    The scansion is determined through linguistic analysis and stored in this function.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    # The scansion pattern derived step-by-step is D S D S D S.
    # D = Dactyl (long-short-short)
    # S = Spondee (long-long)
    scansion = "D S D S D S"

    print(f"Latin Line: {line}")
    print(f"Scansion: {scansion}")

scan_hexameter()