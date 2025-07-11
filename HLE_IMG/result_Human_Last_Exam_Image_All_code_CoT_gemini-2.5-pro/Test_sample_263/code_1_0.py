def analyze_emitter_stability():
    """
    Analyzes the expected stability of three Ir(III) complexes for LEC applications
    by assigning scores based on known stabilizing structural features.
    """
    # --- Step 1: Define Stability Principles ---
    # We assign points for features known to enhance stability:
    # - Steric Hindrance (bulky groups like tert-butyl or large ligands): +2 points
    # - Electronic Stabilization (electron-withdrawing groups like Fluorine): +2 points
    # A baseline ligand gets 0 points.

    # --- Step 2: Evaluate Each Complex ---
    
    # Complex 1: [Ir(ppy)2(bpy)]+
    # - Ancillary ligand (bpy): No special bulky groups.
    # - Cyclometalating ligand (ppy): No electron-withdrawing groups.
    steric_score_1 = 0
    electronic_score_1 = 0
    total_score_1 = steric_score_1 + electronic_score_1
    print("Evaluating Complex 1:")
    print(f"  - Steric hindrance score (from ancillary ligand): {steric_score_1}")
    print(f"  - Electronic stabilization score (from cyclometalating ligand): {electronic_score_1}")
    print(f"  - Final Stability Equation: {steric_score_1} + {electronic_score_1} = {total_score_1}\n")

    # Complex 2: [Ir(ppy)2(large_ligand)]+
    # - Ancillary ligand: Large, bulky, and rigid.
    # - Cyclometalating ligand (ppy): No electron-withdrawing groups.
    steric_score_2 = 2
    electronic_score_2 = 0
    total_score_2 = steric_score_2 + electronic_score_2
    print("Evaluating Complex 2:")
    print(f"  - Steric hindrance score (from ancillary ligand): {steric_score_2}")
    print(f"  - Electronic stabilization score (from cyclometalating ligand): {electronic_score_2}")
    print(f"  - Final Stability Equation: {steric_score_2} + {electronic_score_2} = {total_score_2}\n")

    # Complex 3: [Ir(dfppy)2(dtbbpy)]+
    # - Ancillary ligand (dtbbpy): Has bulky tert-butyl groups.
    # - Cyclometalating ligand (dfppy): Has electron-withdrawing fluorine atoms.
    steric_score_3 = 2
    electronic_score_3 = 2
    total_score_3 = steric_score_3 + electronic_score_3
    print("Evaluating Complex 3:")
    print(f"  - Steric hindrance score (from ancillary ligand): {steric_score_3}")
    print(f"  - Electronic stabilization score (from cyclometalating ligand): {electronic_score_3}")
    print(f"  - Final Stability Equation: {steric_score_3} + {electronic_score_3} = {total_score_3}\n")

    # --- Step 3: Conclude ---
    scores = {
        "Complex 1": total_score_1,
        "Complex 2": total_score_2,
        "Complex 3": total_score_3
    }
    
    most_stable = max(scores, key=scores.get)

    print("--- Conclusion ---")
    print(f"Based on the analysis, {most_stable} achieves the highest stability score.")
    print("It combines both strong steric protection and robust electronic stabilization, which is the most effective strategy for creating stable emitters.")
    print("Therefore, LECs based on Complex 3 are expected to be the most stable.")

if __name__ == '__main__':
    analyze_emitter_stability()