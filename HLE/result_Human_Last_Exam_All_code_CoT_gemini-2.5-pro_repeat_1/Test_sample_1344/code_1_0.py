import textwrap

def propose_catalyst_system():
    """
    Analyzes and proposes an optimal catalyst system for dual-function
    olefin polymerization and polyolefin hydrogenolysis.
    """

    print("--- Catalyst Design Analysis ---")
    print("Objective: A single-site Group IV catalyst for both olefin polymerization and polyolefin hydrogenolysis.")
    print("-" * 34)

    # --- Introduction ---
    intro_text = """
    This task requires a bifunctional catalyst. Olefin polymerization is well-established with Group IV metallocenes. Polyolefin hydrogenolysis (hydrocracking) is more challenging, requiring C-C bond activation. A single-site catalyst must be robust enough for depolymerization conditions while being active for polymerization. The optimal system is likely a supported post-metallocene complex where the support plays a crucial role.
    """
    print("\n[1. Problem Analysis]")
    print(textwrap.dedent(intro_text).strip())

    # --- Metal Selection ---
    metal_text = """
    Group IV metals (Ti, Zr, Hf) are all candidates. Zirconium (Zr) and Hafnium (Hf) are generally preferred over Titanium (Ti) for their higher thermal stability and activity in tandem catalysis. Zirconium offers a good balance of reactivity and cost. We will select Zirconium (Zr) for its proven performance in both polymerization and C-C cleavage reactions.
    """
    print("\n[2. Group IV Metal Selection]")
    print(textwrap.dedent(metal_text).strip())

    # --- Ligand System ---
    ligand_text = """
    The ligand is critical for tuning stability and reactivity. Traditional cyclopentadienyl (Cp) ligands are excellent for polymerization but may lack the robustness for high-temperature hydrogenolysis. A 'pincer' ligand, such as a PNP (Phosphine-Amine-Phosphine) ligand, is an excellent choice. This tridentate ligand offers:
    - High thermal stability due to strong chelation.
    - An 'open' coordination site for reactant binding (H2, olefin).
    - Tunability of electronic properties by modifying the phosphine groups.
    """
    print("\n[3. Ligand System Selection]")
    print(textwrap.dedent(ligand_text).strip())

    # --- Support Selection ---
    support_text = """
    For breaking down a solid polymer, a supported catalyst is essential. The support prevents catalyst aggregation and provides high surface area. A bifunctional approach is optimal, where the support actively participates in the reaction. Sulfated Zirconia (SO4^2-/ZrO2) is a solid superacid support that can:
    - Stabilize the single-site Zr complex.
    - Provide Brønsted acid sites that work with the metal center to heterolytically cleave C-C and C-H bonds in the polymer backbone.
    """
    print("\n[4. Support Selection]")
    print(textwrap.dedent(support_text).strip())

    # --- Proposed Optimal Combination ---
    print("\n[5. Final Proposed Catalyst System]")
    print("Based on the analysis, a highly promising candidate is a Zirconium complex with a PNP pincer ligand, heterogenized on a sulfated zirconia support.")
    print("\n--- Final Equation (System Components) ---")
    print("Metal Center:", "Zr")
    print("Ligand Type:", "PNP Pincer Ligand")
    print("Support:", "Sulfated Zirconia (SO4^2-/ZrO2)")
    print("Mechanism:", "Bifunctional Metal-Acid Catalysis")

if __name__ == '__main__':
    propose_catalyst_system()