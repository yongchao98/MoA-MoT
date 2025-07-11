import itertools

def find_optimal_catalyst():
    """
    This script models the search for an optimal bifunctional catalyst.
    It evaluates combinations of metals, ligands, and supports based on
    heuristic scores for polymerization, hydrogenolysis, and stability.
    """

    # --- 1. Define Components and Heuristic Scores ---
    # Scores are simplified representations (1-10) based on literature trends.
    metals = {
        "Titanium (Ti)": {"polymerization": 6, "hydrogenolysis": 7, "stability": 7},
        "Zirconium (Zr)": {"polymerization": 8, "hydrogenolysis": 8, "stability": 8},
        "Hafnium (Hf)": {"polymerization": 9, "hydrogenolysis": 6, "stability": 9},
    }

    ligands = {
        "Metallocene (e.g., Cp*2)": {"polymerization": 10, "hydrogenolysis": 3, "stability": 8},
        "Constrained-Geometry (CGC)": {"polymerization": 8, "hydrogenolysis": 9, "stability": 7},
        "Phenoxy-Imine (FI)": {"polymerization": 9, "hydrogenolysis": 5, "stability": 6},
        "Bis(phenolate)": {"polymerization": 7, "hydrogenolysis": 6, "stability": 9},
    }

    supports = {
        "None (Homogeneous)": {"polymerization": 1, "hydrogenolysis": 1, "stability": -2},
        "Silica (SiO2)": {"polymerization": 0, "hydrogenolysis": 0, "stability": 2},
        "Alumina (Al2O3)": {"polymerization": -1, "hydrogenolysis": 1, "stability": 1},
    }

    # --- 2. Define Weights for Scoring ---
    # The bifunctional nature requires a balance. Hydrogenolysis is the key challenge.
    weights = {
        "polymerization": 0.3,
        "hydrogenolysis": 0.5,
        "stability": 0.2,
    }

    # --- 3. Iterate and Evaluate All Combinations ---
    best_combo = None
    max_score = -1

    all_components = [metals.items(), ligands.items(), supports.items()]
    # itertools.product generates the Cartesian product of the component lists
    for (metal_name, m_scores), (ligand_name, l_scores), (support_name, s_scores) in itertools.product(*all_components):

        # Calculate the total score for each category
        score_poly = m_scores["polymerization"] + l_scores["polymerization"] + s_scores["polymerization"]
        score_hydro = m_scores["hydrogenolysis"] + l_scores["hydrogenolysis"] + s_scores["hydrogenolysis"]
        score_stab = m_scores["stability"] + l_scores["stability"] + s_scores["stability"]

        # Calculate the final weighted score
        final_score = (weights["polymerization"] * score_poly +
                       weights["hydrogenolysis"] * score_hydro +
                       weights["stability"] * score_stab)

        if final_score > max_score:
            max_score = final_score
            best_combo = {
                "Metal": metal_name,
                "Ligand": ligand_name,
                "Support": support_name,
                "Score": final_score,
                "Score_Poly": score_poly,
                "Score_Hydro": score_hydro,
                "Score_Stab": score_stab,
                "m_scores": m_scores,
                "l_scores": l_scores,
                "s_scores": s_scores
            }

    # --- 4. Report the Best Candidate ---
    if best_combo:
        print("--- Optimal Catalyst Combination based on Model ---")
        print(f"Metal: {best_combo['Metal']}")
        print(f"Ligand: {best_combo['Ligand']}")
        print(f"Support: {best_combo['Support']}")
        print(f"\nCalculated Suitability Score: {best_combo['Score']:.2f}")
        
        # --- 5. Display the "Equation" with numbers ---
        print("\n--- Final Score Calculation ---")
        
        # Polymerization part
        w_p = weights['polymerization']
        m_p = best_combo['m_scores']['polymerization']
        l_p = best_combo['l_scores']['polymerization']
        s_p = best_combo['s_scores']['polymerization']
        print(f"Polymerization Score Contribution: {w_p} * ({m_p} + {l_p} + {s_p}) = {w_p * best_combo['Score_Poly']:.2f}")

        # Hydrogenolysis part
        w_h = weights['hydrogenolysis']
        m_h = best_combo['m_scores']['hydrogenolysis']
        l_h = best_combo['l_scores']['hydrogenolysis']
        s_h = best_combo['s_scores']['hydrogenolysis']
        print(f"Hydrogenolysis Score Contribution: {w_h} * ({m_h} + {l_h} + {s_h}) = {w_h * best_combo['Score_Hydro']:.2f}")

        # Stability part
        w_s = weights['stability']
        m_s = best_combo['m_scores']['stability']
        l_s = best_combo['l_scores']['stability']
        s_s = best_combo['s_scores']['stability']
        print(f"Stability Score Contribution: {w_s} * ({m_s} + {l_s} + {s_s}) = {w_s * best_combo['Score_Stab']:.2f}")

        print("\nFinal Equation:")
        print(f"Total Score = ({w_p} * {best_combo['Score_Poly']}) + ({w_h} * {best_combo['Score_Hydro']}) + ({w_s} * {best_combo['Score_Stab']}) = {best_combo['Score']:.2f}")


if __name__ == "__main__":
    find_optimal_catalyst()