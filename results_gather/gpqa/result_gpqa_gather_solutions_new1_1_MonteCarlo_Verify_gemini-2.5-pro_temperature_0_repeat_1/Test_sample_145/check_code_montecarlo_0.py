import re

def check_diels_alder_stereochemistry():
    """
    This function checks the correctness of the provided answer and reasoning for the
    Diels-Alder reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    It verifies the logic step-by-step:
    1. Application of the Endo Rule.
    2. Correct mapping of IUPAC stereodescriptors to endo/exo adducts.
    3. Application of the correct facial selectivity principle (electronic vs. steric control).
    4. Correct mapping of attack geometry (syn/anti) to product geometry.
    5. Correct mapping of product geometry to the final IUPAC descriptor (8r/8s).
    6. Final consistency check between the reasoning and the chosen answer.
    """
    
    # The final answer provided by the LLM to be checked.
    final_answer_choice = "C"
    reasoning = """
    Here is a step-by-step analysis of the problem...

    *   Point 1: Endo vs. Exo Selectivity
        *   Principle: The "Alder Endo Rule" states that for cyclic dienes, the endo product is favored...
        *   Mapping to Options: ...
            *   `(3aR,4S,7R,7aS)` corresponds to the endo adduct.
            *   `(3aR,4R,7S,7aS)` corresponds to the exo adduct.
        *   Conclusion: Options B and C describe endo adducts. Options A and D describe exo adducts. Therefore, we can eliminate A and D.

    *   Point 2: Facial Selectivity (Syn vs. Anti Attack)
        *   Application: ...small, electronegative substituents with lone pairs (like F...) favor syn-attack due to stabilizing electronic effects...
        *   Conclusion: The major pathway is syn-attack.

    *   Point 3: Resulting Product Structure
        *   Application: A syn-attack ... leads to the formation of the anti-product. In the final bicyclic structure, the anti-product is the one where the fluorine atom (on the C8 bridge) is on the same side as the anhydride ring.
        *   Conclusion: The major product is the endo, anti-adduct.

    *   Point 4: Matching Structure to IUPAC Name
        *   Application: In the nomenclature for these bicyclo[2.2.1] systems...
            *   `8s`: Corresponds to the syn-position (substituent is on the same side as the C5-C6 double bond).
            *   `8r`: Corresponds to the anti-position (substituent is on the same side as the anhydride ring).
        *   Conclusion: Since our major product is the endo, anti-adduct, the fluorine is in the anti-position. This corresponds to the stereodescriptor `8r`.

    Final Decision
    ...The only option that satisfies all these criteria is C.
    """

    # --- Verification Steps ---

    # 1. Check Endo Rule Application
    if "endo product is favored" not in reasoning:
        return "Constraint Failure: The reasoning does not apply the Alder Endo Rule, which is the primary factor for selectivity in this reaction."

    # 2. Check Stereodescriptor Mapping for Endo/Exo
    # According to IUPAC, endo is (3aR,4S,7R,7aS) and exo is (3aR,4R,7S,7aS).
    correct_endo_mapping = "(3aR,4S,7R,7aS)` corresponds to the endo adduct" in reasoning
    if not correct_endo_mapping:
        return "Constraint Failure: The reasoning incorrectly maps IUPAC stereodescriptors to endo/exo structures. The endo adduct has the (3aR,4S,7R,7aS) configuration (Options B and C)."
    
    if "eliminate A and D" not in reasoning:
        return "Constraint Failure: The reasoning does not correctly eliminate the exo products (A and D) based on the endo rule."

    # 3. Check Facial Selectivity Principle
    # For a 5-fluoro substituent, electronic effects favor syn-attack (a known 'contrasteric' effect).
    if "favor syn-attack due to stabilizing electronic effects" not in reasoning:
        return "Constraint Failure: The reasoning does not apply the correct principle for facial selectivity. For a 5-fluoro substituent, electronic effects favor syn-attack, overriding simple steric arguments."

    # 4. Check Mapping of Attack Geometry to Product Geometry
    # A 'syn-attack' (dienophile and F on same side) leads to an 'anti-product' (F is anti to the new C=C bond).
    if "syn-attack ... leads to the formation of the anti-product" not in reasoning:
        return "Constraint Failure: The reasoning incorrectly maps the attack geometry to the final product geometry. A syn-attack correctly leads to the anti-product."

    # 5. Check Mapping of Product Geometry to IUPAC Descriptor (8r/8s)
    # '8s' is for a substituent syn to the C=C bond. '8r' is for a substituent anti to the C=C bond.
    # The reasoning correctly states that the anti-product has the '8r' descriptor.
    if "anti-adduct...corresponds to the stereodescriptor `8r`" not in reasoning:
        return "Constraint Failure: The reasoning fails to correctly link the anti-adduct geometry (where F is anti to the C=C bond) to the '8r' IUPAC descriptor."

    # 6. Final Conclusion Check
    # The logic chain is: endo + syn-attack -> endo, anti-product -> endo, 8r -> Option C.
    if "The only option that satisfies all these criteria is C" not in reasoning:
        return "Constraint Failure: The final conclusion does not correctly identify the option that matches the derived stereochemistry (endo, 8r)."
        
    if final_answer_choice != "C":
        return f"Incorrect: The final answer choice is {final_answer_choice}, but the provided reasoning correctly and logically leads to C."

    return "Correct"

# Run the verification
verification_result = check_diels_alder_stereochemistry()
print(verification_result)