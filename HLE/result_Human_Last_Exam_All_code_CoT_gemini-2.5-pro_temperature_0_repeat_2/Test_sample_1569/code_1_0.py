def diagnose_tick_borne_illness():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    """
    # Step 1: Define the patient's clinical presentation
    patient_data = {
        "age": 27,
        "symptoms": ["fever", "headaches", "myalgia", "disorientation", "heart murmur"],
        "history": "Camping trip to Oklahoma",
        "labs": "Elevated IgM with negative IgG Lyme serology"
    }

    print("--- Patient Case Analysis ---")
    print(f"A {patient_data['age']}-year-old patient presents with: {', '.join(patient_data['symptoms'])}.")
    print(f"Key History: {patient_data['history']}.")
    print(f"Lab Findings: {patient_data['labs']}.")
    print("-" * 29 + "\n")

    # Step 2: Evaluate the differential diagnosis
    print("--- Evaluating Potential Diagnoses ---")

    # A. Babesia microti
    print("A. Babesia microti:")
    print("   - Geography: Primarily Northeast and Upper Midwest US. Less common in Oklahoma.")
    print("   - Conclusion: Poor geographical match.")

    # B. Plasmodium
    print("\nB. Plasmodium (Malaria):")
    print("   - Geography: Not endemic in Oklahoma. Requires travel to specific international regions.")
    print("   - Conclusion: Inconsistent with patient history.")

    # C. Borrelia burgdorferi (Lyme Disease)
    print("\nC. Borrelia burgdorferi (Lyme Disease):")
    print("   - Symptoms: Can cause neurologic and cardiac symptoms.")
    print("   - Labs: The positive IgM/negative IgG fits early Lyme disease.")
    print("   - Geography: Oklahoma is on the fringe of the endemic area.")
    print("   - Conclusion: Possible, but geography is not a strong match.")

    # D. Ehrlichia
    print("\nD. Ehrlichia (Ehrlichiosis):")
    print("   - Geography: Highly endemic in Oklahoma (transmitted by Lone Star tick).")
    print("   - Symptoms: Classic presentation includes fever, headache, myalgia, and confusion.")
    print("   - Labs: Known to cause false-positive Lyme IgM serology.")
    print("   - Conclusion: Excellent match for geography, symptoms, and lab findings.")

    # E. Rickettsia rickettsii (RMSF)
    print("\nE. Rickettsia rickettsii (Rocky Mountain Spotted Fever):")
    print("   - Geography: Highly endemic in Oklahoma.")
    print("   - Symptoms: Similar to Ehrlichiosis, but a characteristic rash is common (not mentioned here).")
    print("   - Conclusion: A strong possibility, but Ehrlichia better explains the specific lab result.")
    print("-" * 36 + "\n")

    # Step 3: Final Conclusion
    print("--- Final Conclusion ---")
    print("The patient's recent travel to Oklahoma, a high-risk area for Ehrlichiosis, is a key clue.")
    print("The clinical syndrome of fever, headache, myalgia, and disorientation is classic for Ehrlichiosis.")
    print("Most importantly, Ehrlichiosis is a well-known cause of a false-positive IgM Lyme test.")
    print("Therefore, the combination of geography, symptoms, and specific lab findings makes Ehrlichia the most likely diagnosis.")

diagnose_tick_borne_illness()
<<<D>>>