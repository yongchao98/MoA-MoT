def find_chromosomal_abnormality():
    """
    This function analyzes a patient's symptoms to determine the most likely
    chromosomal abnormality from a list of options.
    """
    patient_symptoms = [
        "cleft palate",
        "midface hypoplasia",
        "micrognathia",
        "dysplastic ears",
        "prominent eyes",
        "intellectual disability",
        "incomprehensible speech",
        "developmental delay"
    ]

    print("Step 1: Identify the patient's key features.")
    print(f"The patient presents with a distinct constellation of symptoms, including: {', '.join(patient_symptoms)}.")
    print("\nStep 2: Evaluate the options based on the features.")
    print("Equation: Patient's Features vs. Known Syndromes")
    
    print("\n- Analyzing Chromosome 21 (Down Syndrome): Typically associated with a flat facial profile and upslanting palpebral fissures, not the midface hypoplasia and prominent eyes described. This is an unlikely match.")
    
    print("\n- Analyzing Chromosome 13 (Patau Syndrome): This is a severe condition, and survival to age 15 is extremely rare. It is not a likely match.")

    print("\n- Analyzing Chromosome 22 (22q11.2 Deletion Syndrome): This syndrome, also known as Velocardiofacial syndrome, is classically associated with many of the patient's features.")
    print("  - Match 1: Cleft palate is a very common feature.")
    print("  - Match 2: The facial features (midface hypoplasia, micrognathia, prominent eyes, dysplastic ears) are highly characteristic.")
    print("  - Match 3: Intellectual disability and speech/language deficits are core components.")
    print("  - Conclusion: This is a very strong match.")

    print("\nStep 3: Final Conclusion.")
    print("The combination of cleft palate, specific facial dysmorphisms, and developmental/speech delay points strongly to an abnormality on Chromosome 22.")

find_chromosomal_abnormality()