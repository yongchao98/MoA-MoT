def scan_hexameter_line():
    """
    Scans the Latin hexameter line: "verum ubi equi atque hominis casu convenit imago."
    This specific line from Lucretius is known for its complexity, requiring several poetic licenses.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # The feet are determined by a complex analysis involving:
    # 1. Elision: e.g., ver(um) ubi, ub(i) equi, atqu(e) hominis
    # 2. Syncope: hominis is scanned as homnis (--)
    # 3. Tmesis & Verb Tense: convenit is split into 'con' and the perfect 'venit'
    feet_syllables = [
        "verum ubi e-",
        "-qui atque",
        "hominis",
        "casu con-",
        "venit i-",
        "mago"
    ]
    
    scansion_symbols = ["D", "S", "S", "S", "D", "S"]
    
    scansion_explanations = [
        "(Dactyl: — U U)",
        "(Spondee: — —)",
        "(Spondee: — — , via syncope 'homnis')",
        "(Spondee: — — , via tmesis 'con...')",
        "(Dactyl: — U U , via tmesis '...venit')",
        "(Spondee: — —)"
    ]

    print(f"Scansion of the line: \"{line}\"\n")
    
    for i in range(6):
        print(f"Foot {i+1}: {feet_syllables[i]:<15} {scansion_symbols[i]} {scansion_explanations[i]}")
        
    print("\nFinal scansion pattern:")
    
    # Print each symbol of the final pattern as requested
    final_pattern = " ".join(scansion_symbols)
    print(final_pattern)

scan_hexameter_line()