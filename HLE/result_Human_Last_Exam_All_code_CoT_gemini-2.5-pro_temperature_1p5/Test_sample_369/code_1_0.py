def print_scansion():
    """
    Prints the scansion of the Latin hexameter line:
    "verum ubi equi atque hominis casu convenit imago."
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # Based on metrical analysis, the line is divided into feet as follows.
    # Elisions are noted: ver(um) ubi, equ(i) atque, convenit i(mago)
    # The scansion requires some standard poetic licenses.
    
    foot1_word = "verum ubi"
    foot1_syllables = "vē - rŭ - bĭ"
    foot1_meter = "D"
    
    foot2_word = "equi atque"
    foot2_syllables = "ē - qu'ātque"
    foot2_meter = "S"
    
    foot3_word = "hominis"
    foot3_syllables = "hōm - nĭs"
    foot3_meter = "S"

    foot4_word = "casu"
    foot4_syllables = "cā - sū"
    foot4_meter = "S"
    
    foot5_word = "convenit i-"
    foot5_syllables = "cŏn - vĕ - nĭt" #This foot presents a known difficulty.
                                   # For a more regular scansion let's assume convēnit (perfect)
    foot5_syllables_alt = "vē - nĭt - ĭ"
    foot5_meter_alt = "D"
    
    foot6_word = "-mago"
    foot6_syllables = "mā - gō"
    foot6_meter = "S"
    
    # A highly plausible scansion is D-S-S-S-D-S. Let's represent that.
    
    final_feet = {
        1: ("verum ubi", "D"),
        2: ("equi atque", "S"),
        3: ("hominis", "S"),
        4: ("casu", "S"),
        5: ("convenit i-", "D"),
        6: ("-mago", "S")
    }

    print(f"Original line: {line}\n")
    print("Scansion:")
    
    # This reconstruction shows a plausible division of the line.
    # It requires assuming ver(um) ubi elision, equ(i) atque elision (with a poetically lengthened 'e'),
    # hominis undergoes syncope (hom'nis), and convenit i- forms a dactyl with the 'i' of imago.
    
    scanned_line = "vērum‿ub(i)‿ēqu(i)‿ātque | hōmnis | cāsū | cōnvĕnĭt‿ĭ | māgō"
    scansion_pattern = "D S S D D S" 
    
    # Re-evaluating the most standard pattern. Let's use D S S S D S.
    
    final_pattern = "D S S S D S"
    final_division = "verum‿ubi | equi‿atque | hominis | cāsū | convēnit‿i|māgō"
    
    print(final_division)
    print(final_pattern)


print_scansion()
<<<D S S S D S>>>