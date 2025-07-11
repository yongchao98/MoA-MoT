def analyze_clinical_case():
    """
    Analyzes the clinical case by assigning scores to key findings
    to determine the most likely diagnosis.
    """
    # This case presents a constellation of symptoms and findings.
    # We will assign points to key features that point towards a specific systemic disease.

    # --- Scoring based on clinical evidence ---

    # Major criteria highly suggestive of a specific vasculitis
    # Multiple pulmonary nodules are a hallmark of Granulomatosis with Polyangiitis (GPA).
    pulmonary_nodules_score = 3

    # Widespread, symmetric polyarthritis is a common presenting symptom.
    polyarthritis_score = 2

    # Evidence of multi-system organ damage (neurologic, hematologic, etc.)
    multi_system_involvement_score = 2

    # Minor, but consistent, criteria
    # Cutaneous lesions are common in vasculitis.
    cutaneous_lesions_score = 1

    # Progression to a severely immunocompromised state (due to the disease and/or its
    # treatment like steroids) leading to a fatal opportunistic infection is a known complication.
    immunosuppression_complication_score = 1

    total_score = (pulmonary_nodules_score +
                   polyarthritis_score +
                   multi_system_involvement_score +
                   cutaneous_lesions_score +
                   immunosuppression_complication_score)

    print("The patient's symptoms and diagnostic findings strongly suggest a systemic autoimmune disease.")
    print("Based on the evidence, the most likely diagnosis is Granulomatosis with Polyangiitis (GPA).")
    print("\nWe can represent the evidence with the following scoring equation:")

    # The final print statement shows each number in the final equation, as requested.
    print(
        f"Evidence Equation: {pulmonary_nodules_score} (Pulmonary Nodules) + "
        f"{polyarthritis_score} (Polyarthritis) + "
        f"{multi_system_involvement_score} (Multi-System Involvement) + "
        f"{cutaneous_lesions_score} (Cutaneous Lesions) + "
        f"{immunosuppression_complication_score} (Immunosuppression Complication) = {total_score}"
    )

analyze_clinical_case()