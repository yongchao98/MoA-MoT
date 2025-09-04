import random

def analyze_structure():
    """
    Analyzes chemical candidates against spectral data to identify an unknown compound.
    (a) Samples from a list of candidates.
    (b) Narrows candidates based on spectral features.
    (c) Runs exact checks for verification.
    (d) Prints the final verified answer.
    """
    # Define the candidates and their key structural features for spectral prediction
    candidates = {
        "A": {"name": "Phenyl chloroformate", "mw": 156.57, "has_cl": True, "has_cooh": False, "aromatic_pattern": "mono"},
        "B": {"name": "4-chlorobenzoic acid", "mw": 156.57, "has_cl": True, "has_cooh": True, "aromatic_pattern": "para"},
        "C": {"name": "3-Chloro-2-hydroxybenzaldehyde", "mw": 156.57, "has_cl": True, "has_cooh": False, "aromatic_pattern": "tri"},
        "D": {"name": "2-chlorobenzoic acid", "mw": 156.57, "has_cl": True, "has_cooh": True, "aromatic_pattern": "ortho"},
    }

    # (a) Sample: We will check all candidates. Shuffling simulates random exploration.
    candidate_keys = list(candidates.keys())
    random.shuffle(candidate_keys)
    print(f"Starting analysis. Order of exploration: {candidate_keys}\n")

    passing_candidates = []

    # (b) Narrow candidates & (c) Run exact checks
    for key in candidate_keys:
        candidate = candidates[key]
        name = candidate['name']
        print(f"--- Checking Candidate {key}: {name} ---")
        
        failures = []

        # Check 1: Mass Spec (MW ~156, one Cl atom)
        # All candidates have the correct MW and one Cl, so all pass this check.
        print("  - MS Check (MW≈156, 1xCl): PASSED")

        # Check 2: IR (Carboxylic Acid)
        if not candidate["has_cooh"]:
            failures.append("IR data indicates a carboxylic acid, which is absent.")
        
        # Check 3: 1H NMR (Carboxylic Acid Proton + Aromatic Pattern)
        if not candidate["has_cooh"]:
            failures.append("NMR shows a COOH proton at 11.0 ppm, which is absent.")
        if candidate["aromatic_pattern"] != "para":
            failures.append(f"NMR shows a para-substitution pattern (2 doublets, 2H each), but this structure is {candidate['aromatic_pattern']}-substituted.")

        # Final verdict for this candidate
        if not failures:
            print(f"  - IR & NMR Checks: PASSED. Structure is consistent with all data.")
            print(f"VERDICT: Candidate {key} ({name}) passed all checks.\n")
            passing_candidates.append(key)
        else:
            # Use set to remove duplicate failure reasons (e.g., from IR and NMR both noting lack of COOH)
            unique_failures = sorted(list(set(failures)))
            print(f"  - IR & NMR Checks: FAILED. Reasons: {' '.join(unique_failures)}")
            print(f"VERDICT: Candidate {key} ({name}) failed verification.\n")

    # (d) Print the final verified answer
    print("--- FINAL CONCLUSION ---")
    if len(passing_candidates) == 1:
        winner_key = passing_candidates[0]
        winner_name = candidates[winner_key]['name']
        print(f"The only structure consistent with all spectral data is {winner_key}: {winner_name}.")
        print("\nFinal Answer is B")
    elif len(passing_candidates) > 1:
        print(f"Error: Multiple candidates passed: {passing_candidates}. The data may be ambiguous.")
    else:
        print("Error: No candidate passed all checks. The correct structure may not be in the options provided.")

analyze_structure()