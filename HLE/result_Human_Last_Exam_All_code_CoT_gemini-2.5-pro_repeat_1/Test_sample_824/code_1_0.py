def identify_key_lab_parameter():
    """
    Analyzes the clinical scenario to determine the most indicative lab parameter
    for the cause of rapid renal decline.
    """

    # Patient's history and symptoms point towards Systemic Lupus Erythematosus (SLE).
    # Key features: facial rash, joint pain, fever, kidney involvement.
    patient_diagnosis = "Systemic Lupus Erythematosus (SLE) with Lupus Nephritis"

    # The cause of kidney damage in Lupus Nephritis is a specific immunological process.
    cause_of_damage = "Deposition of immune complexes in the kidneys"

    # This process triggers a cascade that consumes specific immune proteins.
    immunological_consequence = "Activation and consumption of the complement system"

    # The measurable result of this consumption is a key lab finding.
    key_lab_indicator = "Low serum complement levels (C3 and C4)"

    print("Clinical Analysis:")
    print(f"The patient's symptoms strongly suggest a severe flare of: {patient_diagnosis}")
    print(f"The direct cause of the kidney damage in this condition is: {cause_of_damage}")
    print(f"This leads to the following immunological event: {immunological_consequence}")
    print("\nTherefore, the lab parameter that best indicates the cause of the rapid renal function decline is a finding that reflects this event.")

    # The prompt requests an 'equation' format. We will represent the diagnostic logic this way.
    print("\nDiagnostic Logic Equation:")
    part_1 = "Active SLE Flare"
    part_2 = "Immune Complex Deposition"
    part_3 = "Complement Consumption"
    final_answer = key_lab_indicator

    # Printing each "number" or component of the final "equation"
    print(f"Step 1: {part_1}")
    print(f"Step 2: {part_2}")
    print(f"Step 3: {part_3}")
    print(f"Conclusion: The best indicator is -> {final_answer}")


if __name__ == "__main__":
    identify_key_lab_parameter()
