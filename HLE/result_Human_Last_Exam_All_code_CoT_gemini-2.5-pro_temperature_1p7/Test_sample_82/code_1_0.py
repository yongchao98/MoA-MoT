import sys

def analyze_heritability():
    """
    Analyzes the statements about heritability and polygenic scores
    based on the provided information.
    """
    # --- Part 1: Define parameters from the problem ---
    # We are given broad-sense heritability (H²)
    H2 = 0.5

    # For demonstration, let's assume a total phenotypic variance (Vp)
    Vp = 100.0

    # From H², we can calculate the total genetic variance (Vg)
    Vg = H2 * Vp

    print("--- Foundational Equations and Known Values ---")
    print(f"Broad-sense heritability, H² = Vg / Vp = {H2}")
    print(f"Assuming Vp = {Vp}, total genetic variance Vg = {H2} * {Vp} = {Vg}")
    print("Genetic variance, Vg, is composed of Va (additive), Vd (dominance), and Vi (epistatic).")
    print("A standard linear Polygenic Score (PGS) can, at best, explain the variance due to Va.")
    print("The maximum variance explained by a PGS is the narrow-sense heritability, h² = Va / Vp.\n")

    # --- Part 2: Analysis of Statement A ---
    print("--- Analysis of Statement A ---")
    print("A. The polygenic score can not explain more than 50% of the variance.")
    print("By definition, Va <= Vg.")
    print(f"Dividing by Vp gives: Va / Vp <= Vg / Vp.")
    print(f"This translates to: h² <= H².")
    print(f"Since H² = {H2}, the maximum possible value for h² is {H2}.")
    print("The variance explained by the PGS cannot exceed 50%.")
    print("Conclusion: Statement A is necessarily true.\n")

    # --- Part 3: Analysis of Statement C ---
    print("--- Analysis of Statement C ---")
    print("C. ...the PGS... will not approach... 50% due to... non-linear effects...")
    print("This statement implies that h² will be strictly less than H².")
    print("This happens if there is any non-additive genetic variance (Vd > 0 or Vi > 0).")
    print("Let's model a realistic scenario for a polygenic trait:")
    
    # Assume some Vg is non-additive, e.g., 80% is additive, 20% is non-additive
    Va = Vg * 0.8
    V_non_additive = Vg * 0.2
    
    # Calculate the resulting h²
    h2 = Va / Vp

    print(f"  Assume Va = {Va:.1f} and non-additive variance (Vd + Vi) = {V_non_additive:.1f}.")
    print(f"  This is a realistic scenario where Va ({Va:.1f}) + V_non_additive ({V_non_additive:.1f}) = Vg ({Vg:.1f}).")
    print(f"  In this case, the final equation for h² is:")
    print(f"  h² = Va / Vp = {Va:.1f} / {Vp:.1f} = {h2}")
    print(f"The result h² = {h2} is less than H² = {H2}.")
    print("Since 'polygenic' implies multiple genes, the existence of such non-linear interactions is a standard biological assumption.")
    print("Conclusion: Statement C is considered true.\n")

    # --- Part 4: Final Answer ---
    print("--- Overall Conclusion ---")
    print("Both statements A and C are correct.")
    print("A is true by mathematical definition, and C is true under standard biological assumptions for polygenic traits.")
    print("Therefore, the correct choice is E.")

if __name__ == "__main__":
    analyze_heritability()
<<<E>>>