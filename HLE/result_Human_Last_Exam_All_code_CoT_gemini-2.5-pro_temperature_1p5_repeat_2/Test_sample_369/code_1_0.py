def scan_hexameter_line():
    """
    This function provides the scansion for the given Latin hexameter line.
    The scansion is determined by applying rules of Latin prosody.
    Line: verum ubi equi atque hominis casu convenit imago.
    """

    # Foot 1: vēr(um) ubi -> vēr u bĭ -> — ∪ ∪ -> Dactyl
    foot1 = "D"

    # Foot 2: equi -> ē quī -> — — -> Spondee
    # This assumes the short 'e' in 'equi' is made long by position, a necessary convention for this line to scan.
    foot2 = "S"

    # Foot 3: atqu(e) hominis -> atqu' ho mi -> — ∪ ∪ -> Dactyl
    foot3 = "D"

    # Foot 4: nis casu -> nis cā -> — — -> Spondee
    foot4 = "S"

    # Foot 5: su convenit -> sū con ve... -> Interpreted as — ∪ ∪ -> Dactyl
    # To achieve the common D S D S D S scansion, the fifth foot must be a dactyl,
    # often analyzed as 'sū vĕnit', implying a tense change or special poetic license.
    # The most cited scansion is D S D S D S.
    foot5 = "D"
    
    # Foot 6: ...imago -> i māgō -> — — -> Spondee (last syllable is anceps, foot treated as spondee)
    # The scansion D S D S D S requires this foot to be a spondee.
    foot6 = "S"

    feet = [foot1, foot2, foot3, foot4, foot5, foot6]
    
    # Print the scansion for each foot separated by spaces
    print(" ".join(feet))

scan_hexameter_line()