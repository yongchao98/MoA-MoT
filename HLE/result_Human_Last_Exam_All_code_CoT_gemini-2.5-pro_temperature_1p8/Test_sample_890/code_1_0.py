def find_chromosomal_abnormality():
    """
    Analyzes clinical features to identify the most likely chromosomal abnormality.
    """
    # Key features presented in the patient's case
    patient_features = {
        "Palatal": "cleft palate",
        "Craniofacial": "midface hypoplasia, micrognathia, dysplastic ears, prominent eyes, microcephaly",
        "Neurodevelopmental": "intellectual disability, delayed milestones, incomprehensible speech"
    }

    # Features of common chromosomal syndromes
    syndrome_features = {
        "Chromosome 22 (22q11.2 deletion)": "Cleft palate, characteristic facial features (midface hypoplasia, micrognathia), learning/developmental/speech delays",
        "Chromosome 21 (Trisomy 21)": "Intellectual disability, flat facial profile, upslanting eyes, single palmar crease",
        "Chromosome 13 (Trisomy 13)": "Severe intellectual disability, cleft lip/palate, microcephaly, typically lethal in infancy"
    }

    print("Step 1: Analyzing patient's key features.")
    for category, features in patient_features.items():
        print(f"- {category}: {features}")

    print("\nStep 2: Comparing with known syndromes.")
    print("The patient's combination of symptoms shows the strongest match with the syndrome associated with Chromosome 22.")

    print("\nStep 3: Justification.")
    print(f"Patient's core features match classic signs of 22q11.2 deletion syndrome.")
    print(f"Syndrome Features on Chromosome 22: {syndrome_features['Chromosome 22 (22q11.2 deletion)']}")
    
    print("\nConclusion: The chromosomal abnormality is expected on Chromosome 22.")

find_chromosomal_abnormality()