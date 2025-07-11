def solve_entailment():
    """
    Calculates the final projected natural logic operator for the given inference
    by simulating the composition of semantic relations according to MacCartney's framework.
    """
    # Define the 7 semantic relations and their symbols for printing
    RELATIONS = {
        "EQU": ("Equivalence", "≡"),
        "FE": ("Forward Entailment", "⊏"),
        "RE": ("Reverse Entailment", "⊐"),
        "NEG": ("Negation", "^"),
        "ALT": ("Alternation", "|"),
        "COV": ("Cover", "v"),
        "IND": ("Independence", "#"),
    }
    
    # Define the order of relations for the table index
    REL_ORDER = ["EQU", "FE", "RE", "NEG", "ALT", "COV", "IND"]
    
    # MacCartney's composition table C(row, col)
    COMPOSITION_TABLE = [
        # EQU   FE     RE     NEG    ALT    COV    IND
        ["EQU", "FE",  "RE",  "NEG", "ALT", "COV", "IND"], # EQU
        ["FE",  "FE",  "IND", "COV", "ALT", "COV", "IND"], # FE
        ["RE",  "IND", "RE",  "ALT", "ALT", "IND", "IND"], # RE
        ["NEG", "RE",  "FE",  "EQU", "COV", "ALT", "IND"], # NEG
        ["ALT", "IND", "ALT", "COV", "ALT", "FE",  "IND"], # ALT
        ["COV", "COV", "IND", "FE",  "COV", "IND", "IND"], # COV
        ["IND", "IND", "IND", "IND", "IND", "IND", "IND"], # IND
    ]

    def compose(rel1_key, rel2_key):
        """Composes two relations using the composition table."""
        idx1 = REL_ORDER.index(rel1_key)
        idx2 = REL_ORDER.index(rel2_key)
        return COMPOSITION_TABLE[idx1][idx2]

    # --- Start of Inference Calculation ---

    print("Step 1: Determine relation R(X, Y) between the core arguments.")
    print("  - Premise Argument (X): 'a pop song by Taylor Swift'")
    print("  - Hypothesis Argument (Y): 'a song by Michael Jackson'")
    print("-" * 20)

    # Edit 1: 'pop song' -> 'song'. This is a generalization.
    # The relation from a specific to a general is Forward Entailment.
    rel_edit1_key = "FE"
    rel_edit1_name, rel_edit1_sym = RELATIONS[rel_edit1_key]
    print(f"Edit 1 (Generalization): 'pop song' -> 'song'")
    print(f"  - Relation: {rel_edit1_name} ({rel_edit1_sym})")

    # The current relation between arguments starts with this first edit.
    current_rel_xy_key = rel_edit1_key

    # Edit 2: 'Taylor Swift' -> 'Michael Jackson'. These are exclusive.
    # The relation is Negation.
    rel_edit2_key = "NEG"
    rel_edit2_name, rel_edit2_sym = RELATIONS[rel_edit2_key]
    print(f"\nEdit 2 (Substitution): 'Taylor Swift' -> 'Michael Jackson'")
    print(f"  - Relation: {rel_edit2_name} ({rel_edit2_sym})")
    
    print(f"\nComposition for R(X, Y):")
    current_rel_xy_name, current_rel_xy_sym = RELATIONS[current_rel_xy_key]
    print(f"  [{current_rel_xy_name} ({current_rel_xy_sym})] ○ [{rel_edit2_name} ({rel_edit2_sym})]")
    
    # Perform composition
    rel_xy_key = compose(current_rel_xy_key, rel_edit2_key)
    rel_xy_name, rel_xy_sym = RELATIONS[rel_xy_key]
    
    print(f"  = {rel_xy_name} ({rel_xy_sym})")
    print("-" * 20)
    
    print("\nStep 2: Apply the final negation from the hypothesis ('is not...').")
    print("  - The relation from premise ('singing X') to the un-negated hypothesis ('singing Y') is R(X, Y).")
    print(f"  - We compose R(X, Y) with the final Negation operator.")

    final_neg_key = "NEG"
    final_neg_name, final_neg_sym = RELATIONS[final_neg_key]

    print(f"\nFinal Composition:")
    print(f"  [{rel_xy_name} ({rel_xy_sym})] ○ [{final_neg_name} ({final_neg_sym})]")
    
    # Final composition
    final_rel_key = compose(rel_xy_key, final_neg_key)
    final_rel_name, final_rel_sym = RELATIONS[final_rel_key]
    
    print(f"  = {final_rel_name} ({final_rel_sym})")
    print("-" * 20)

    print(f"\nThe name of the final projected natural logic operator is: {final_rel_name}")

solve_entailment()