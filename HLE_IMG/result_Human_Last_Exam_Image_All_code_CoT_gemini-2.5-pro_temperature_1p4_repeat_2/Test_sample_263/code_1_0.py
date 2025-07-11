def analyze_complex_stability():
    """
    Analyzes the expected stability of three Ir(III) complexes for LEC applications
    based on their molecular structures.
    """
    # Molecular features known to enhance stability in Ir(III) emitters:
    # 1. Fluorination of the cyclometalating (C^N) ligand: Increases oxidative stability.
    # 2. Bulky groups on the ancillary (N^N) ligand: Provide steric protection, preventing degradation.
    # 3. Rigid N^N ligand: Can offer minor stability improvements over flexible ones.

    complexes = {
        "Complex 1": {
            "C^N Ligand": "ppy (phenylpyridine)",
            "N^N Ligand": "bpy (bipyridine)",
            "Features": []
        },
        "Complex 2": {
            "C^N Ligand": "ppy (phenylpyridine)",
            "N^N Ligand": "phen-imidazole derivative",
            "Features": ["Rigid N^N ligand"]
        },
        "Complex 3": {
            "C^N Ligand": "dfppy (difluorophenylpyridine)",
            "N^N Ligand": "dtbbpy (di-tert-butyl-bipyridine)",
            "Features": ["Fluorinated C^N ligand", "Bulky N^N ligand"]
        }
    }

    print("Step-by-step analysis of complex stability for LEC devices:\n")

    stability_scores = {}
    for name, data in complexes.items():
        score = 0
        print(f"--- Analyzing {name} ---")
        print(f"  C^N Ligand: {data['C^N Ligand']}")
        print(f"  N^N Ligand: {data['N^N Ligand']}")

        if not data["Features"]:
            print("  - This complex is a baseline structure with no special stabilizing modifications.")
            score = 1
        else:
            if "Fluorinated C^N ligand" in data["Features"]:
                print("  - Feature: Fluorination. This significantly increases resistance to electrochemical degradation.")
                score += 2
            if "Bulky N^N ligand" in data["Features"]:
                print("  - Feature: Bulky tert-butyl groups. These provide steric protection, shielding the complex from degradation pathways.")
                score += 2
            if "Rigid N^N ligand" in data["Features"]:
                print("  - Feature: Rigid N^N ligand backbone. This offers a minor improvement in stability compared to simple bipyridine.")
                score += 1
        stability_scores[name] = score
        print(f"  Assigned Stability Score: {score}\n")

    # Find the complex with the highest stability score
    most_stable_complex = max(stability_scores, key=stability_scores.get)
    
    print("--- Conclusion ---")
    print("To achieve high stability in LECs, emitter molecules are often modified with specific functional groups.")
    print("Key strategies include:")
    print("1. Fluorination: To increase electrochemical stability (seen in Complex 3).")
    print("2. Bulky Substituents: To provide steric shielding and prevent degradation (seen in Complex 3).")
    print(f"\nComplex 3 incorporates both of these powerful stabilizing strategies, making it the most robust design.")
    print(f"Therefore, LECs based on {most_stable_complex} are expected to be the most stable.")

analyze_complex_stability()