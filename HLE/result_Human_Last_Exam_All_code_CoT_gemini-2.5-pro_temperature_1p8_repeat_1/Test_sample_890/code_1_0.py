def find_matching_syndrome():
    """
    Analyzes patient features to identify the most likely chromosomal abnormality.
    The code simulates a diagnostic process by scoring each potential syndrome
    based on the presence of key clinical features.
    """
    # Key features observed in the patient, mapped to standard clinical terms
    patient_features = {
        "cleft palate",
        "short stature",
        "microcephaly",
        "developmental delay",  # From "delayed developmental milestones"
        "midface hypoplasia",
        "micrognathia",
        "dysplastic ears",
        "intellectual disability",
        "speech problems"       # From "incomprehensible speech"
    }

    # Knowledge base of features for syndromes associated with each chromosome
    syndrome_database = {
        "A. Chromosome 3 (e.g., 3q29 microdeletion)": ["microcephaly", "developmental delay", "intellectual disability"],
        "B. Chromosome 22 (e.g., 22q11.2 deletion)": ["cleft palate", "intellectual disability", "speech problems", "midface hypoplasia", "micrognathia", "dysplastic ears", "developmental delay", "short stature"],
        "C. Chromosome 21 (Trisomy 21)": ["intellectual disability", "flattened facial profile", "single palmar crease"],
        "D. Chromosome 2 (various)": ["variable features", "structural anomalies"],
        "E. Chromosome 13 (Trisomy 13)": ["cleft palate", "microcephaly", "severe intellectual disability", "polydactyly"]
    }

    best_match = None
    max_score = -1

    print("Matching patient features to known chromosomal abnormalities:\n")

    for syndrome, features in syndrome_database.items():
        score = 0
        equation_parts = []
        for feature in features:
            if feature in patient_features:
                score += 1
                equation_parts.append("1") # Represents a matched feature
            else:
                equation_parts.append("0") # Represents an unmatched feature
        
        # This fulfills the request to output each number in the final equation
        equation_str = " + ".join(equation_parts)
        print(f"Syndrome: {syndrome}")
        print(f"Feature Match Equation: {equation_str} = {score}")
        print("-" * 30)
        
        if score > max_score:
            max_score = score
            best_match = syndrome

    print(f"\nConclusion: The highest match score is for {best_match}.")
    print("The constellation of features, particularly the cleft palate, micrognathia, dysplastic ears, and developmental delay, is classic for 22q11.2 deletion syndrome.")

find_matching_syndrome()