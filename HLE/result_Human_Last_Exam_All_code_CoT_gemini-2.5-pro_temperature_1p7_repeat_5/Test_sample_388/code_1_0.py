def generate_counseling_recommendation():
    """
    Analyzes a patient's medications to identify potential interactions and
    generates a counseling recommendation.
    """
    # Medications Allison is taking
    # Excedrin's active ingredients include Aspirin (an NSAID), Acetaminophen, and Caffeine.
    # The key interaction is between the NSAID and the SSRI.
    prescriptions = {
        "Atorvastatin 20mg": "Statin",
        "Junel Fe 1.5/30mg": "Oral Contraceptive",
        "Fluoxetine 20mg": "SSRI"
    }
    otc_medication = {
        "Excedrin": "Pain Reliever (contains NSAID)"
    }

    all_meds = {**prescriptions, **otc_medication}

    # Check for the specific interaction between an SSRI and an NSAID
    has_ssri = "SSRI" in all_meds.values()
    has_nsaid = "Pain Reliever (contains NSAID)" in all_meds.values()

    print("Analyzing Allison's Medications...")
    for med, m_type in all_meds.items():
        print(f"- {med} ({m_type})")
    print("-" * 40)

    if has_ssri and has_nsaid:
        print("Pharmacist Counseling Recommendation:")
        print("\n[POTENTIAL DRUG INTERACTION IDENTIFIED]")
        print("There is a significant interaction risk between Fluoxetine and the Aspirin in Excedrin.")
        print("\n[THE RISK]")
        print("Fluoxetine, which is an SSRI, and Aspirin, which is an NSAID, both independently increase the risk of gastrointestinal (GI) bleeding. When taken together, this risk is significantly increased.")
        print("\n[RECOMMENDATION]")
        print("1. Advise Allison to avoid using NSAIDs like Aspirin (in Excedrin) or Ibuprofen while she is taking Fluoxetine.")
        print("2. Recommend a safer alternative for her headache, such as a medication containing only Acetaminophen (e.g., Tylenol), as long as she has no other contraindications.")
        print("3. Encourage her to discuss headache management with her doctor to ensure her treatment is safe.")
    else:
        print("No significant interactions found based on the provided information.")

generate_counseling_recommendation()