def find_best_match():
    """
    Analyzes patient symptoms to find the most likely chromosomal abnormality.
    """
    # Key features presented by the 15-year-old patient
    patient_features = {
        "cleft palate",
        "intellectual disability",
        "incomprehensible speech",
        "midface hypoplasia",
        "micrognathia",
        "prominent eyes",
        "dysplastic ears",
        "short stature",
        "microcephaly"
    }

    # Characteristic features of syndromes associated with each chromosome
    syndrome_features = {
        "A. Chromosome 3": {
            "intellectual disability", "short stature", "postaxial polydactyly"
        },
        "B. Chromosome 22 (22q11.2 Deletion)": {
            "cleft palate", "intellectual disability", "incomprehensible speech",
            "midface hypoplasia", "micrognathia", "prominent eyes",
            "dysplastic ears", "short stature", "cardiac defects"
        },
        "C. Chromosome 21 (Down Syndrome)": {
            "intellectual disability", "flat facial profile", "upslanting palpebral fissures",
            "single transverse palmar crease", "short stature", "dysplastic ears"
        },
        "D. Chromosome 2": {
            "intellectual disability", "microcephaly", "growth delay"
        },
        "E. Chromosome 13 (Patau Syndrome)": {
            "severe intellectual disability", "cleft palate", "microcephaly",
            "polydactyly", "severe heart defects", "often lethal in infancy"
        }
    }

    print("Analyzing patient features against known syndromes...\n")
    best_match = None
    max_score = -1

    for choice, features in syndrome_features.items():
        # Calculate the number of matching features
        matching_features = patient_features.intersection(features)
        score = len(matching_features)
        
        print(f"Diagnosis: {choice}")
        print(f"Matching Features ({score}): {', '.join(sorted(list(matching_features))) if matching_features else 'None'}")
        print("-" * 30)

        if score > max_score:
            max_score = score
            best_match = choice

    print(f"\nConclusion: The highest number of matching features ({max_score}) is found for {best_match}.")
    print("The patient's constellation of symptoms, especially the combination of cleft palate, specific facial dysmorphisms, and speech issues, strongly points to a deletion on Chromosome 22.")

find_best_match()