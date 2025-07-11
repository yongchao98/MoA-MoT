def solve_clinical_case():
    """
    This function analyzes the clinical scenario and identifies the key lab parameter.
    
    The patient's history of facial rash, joint pain, and kidney issues points to
    Systemic Lupus Erythematosus (SLE). The rapid decline in kidney function
    after stopping steroids suggests a severe flare of Lupus Nephritis.

    The mechanism of kidney damage in lupus nephritis is the deposition of
    immune complexes in the glomeruli. These complexes activate a part of the
    immune system called the complement cascade.

    During a severe flare, this activation is so intense that it consumes the
    complement proteins, leading to measurably low levels in the blood.
    """
    
    # Key parameters in the complement system that are consumed during a flare.
    protein_1 = 3 # Represents C3
    protein_2 = 4 # Represents C4

    # Explanation of the final answer
    explanation = (
        "The patient's presentation suggests a severe flare of Lupus Nephritis, a complication of Systemic Lupus Erythematosus (SLE).\n\n"
        "In this condition, immune complexes (antibody-antigen clusters) deposit in the kidneys, triggering an intense inflammatory response. "
        "This process activates and consumes proteins of the complement system.\n\n"
        "Therefore, the lab parameter that could have best indicated the cause of the rapid renal function decline is:\n"
        f"Low C{protein_1} and C{protein_2} complement levels.\n\n"
        "A sharp drop in these levels is a direct marker of the massive immune complex activity responsible for destroying the kidney tissue."
    )
    
    print(explanation)

solve_clinical_case()