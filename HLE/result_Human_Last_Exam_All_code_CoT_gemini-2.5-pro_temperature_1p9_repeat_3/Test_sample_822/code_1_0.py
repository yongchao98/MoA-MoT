def diagnose_patient():
    """
    Analyzes the provided medical case study to determine the most likely disease.
    This function structures the diagnostic reasoning step by step.
    """

    # --- Step 1: Identify Key Information from the Case ---
    risk_factors = {
        "Smoking": "20-pack-year history, a major risk factor for lung cancer.",
        "Occupation": "Ship building, implying significant asbestos exposure, another major risk factor for lung cancer."
    }

    chronic_findings = {
        "Initial Symptoms": "Pain/swelling in wrists, ankles, elbows, suggesting hypertrophic pulmonary osteoarthropathy (HPOA), a paraneoplastic syndrome strongly associated with lung cancer.",
        "Chest X-ray": "Multiple pulmonary nodules, a classic finding in lung cancer."
    }

    acute_illness_findings = {
        "Clinical Presentation": "Fever, productive cough, cutaneous lesions in an immunocompromised host.",
        "Treatment Failure": "Ineffective aminoglycoside therapy, suggesting an atypical pathogen.",
        "Likely Pathogen": "This constellation is classic for Nocardiosis, an opportunistic infection."
    }

    # --- Step 2: Print the Reasoning Process ---
    print("Patient Case Analysis:")
    print("=" * 25)

    print("\n[Primary Disease Analysis]")
    print(f"1. Risk Factors: The patient's history includes heavy smoking ({risk_factors['Smoking']})")
    print(f"   and work in ship building ({risk_factors['Occupation']}).")
    print("   This combination overwhelmingly increases the risk for a specific malignancy.")
    print("\n2. Clinical Evidence: The initial symptoms ({chronic_findings['Initial Symptoms']})")
    print("   and the chest X-ray findings ({chronic_findings['Chest X-ray']})")
    print("   strongly point towards this malignancy as the underlying cause.")

    print("\n[Terminal Illness Analysis]")
    print("3. Complication: The patient's underlying cancer made him immunocompromised.")
    print(f"   He developed an opportunistic infection ({acute_illness_findings['Clinical Presentation']})")
    print(f"   that was refractory to standard treatment ({acute_illness_findings['Treatment Failure']}).")
    print(f"   The specific symptoms point to {acute_illness_findings['Likely Pathogen']}.")

    print("\n" + "=" * 25)

    # --- Step 3: Conclude with the Final Diagnosis ---
    final_diagnosis = "Lung Cancer"
    print(f"\nConclusion: The primary disease that explains the entire clinical course, from the initial joint pains to the final fatal infection, is Lung Cancer.")
    print(f"\nThe Nocardiosis was a severe, terminal complication of the cancer.")
    print(f"\nFinal Answer: {final_diagnosis}")


if __name__ == "__main__":
    diagnose_patient()