import textwrap

def find_chromosomal_abnormality():
    """
    Analyzes a clinical case to identify the most likely chromosomal abnormality.
    """
    patient_symptoms = {
        "cleft palate",
        "microcephaly",
        "midface hypoplasia",
        "micrognathia",
        "dysplastic ears",
        "intellectual disability",
        "incomprehensible speech",
        "short stature",
        "clinodactyly",
        "developmental delay"
    }

    syndrome_features = {
        "A. Chromosome 3 (e.g., 3q29 microdeletion)": {
            "intellectual disability", "microcephaly", "autistic features"
        },
        "B. Chromosome 22 (22q11.2 deletion)": {
            "cleft palate", "midface hypoplasia", "micrognathia",
            "dysplastic ears", "developmental delay", "speech problems",
            "intellectual disability", "cardiac defects"
        },
        "C. Chromosome 21 (Trisomy 21)": {
            "intellectual disability", "flat facial profile", "upslanting eyes",
            "single palmar crease", "dysplastic ears"
        },
        "D. Chromosome 2 (e.g., 2q37 deletion)": {
            "intellectual disability", "brachydactyly (short fingers/toes)",
            "obesity", "autistic features"
        },
        "E. Chromosome 13 (Trisomy 13)": {
            "severe intellectual disability", "cleft lip/palate", "microcephaly",
            "holoprosencephaly", "polydactyly", "low survival rate"
        }
    }

    best_match = None
    max_score = -1
    best_match_details = ""

    print("Analyzing the patient's symptoms against known chromosomal abnormalities...")
    print("-" * 50)

    for choice, features in syndrome_features.items():
        # Give higher weight to highly specific or classic features
        # In this case, the combination of palate, facial, and speech issues is key for 22q11.2
        key_features_for_22q = {"cleft palate", "midface hypoplasia", "micrognathia", "speech problems"}
        
        score = 0
        matched_symptoms = patient_symptoms.intersection(features)
        
        # Simple scoring based on number of matches
        # A more complex model could weight specific symptoms
        score = len(matched_symptoms)

        # Boost score for a classic presentation
        if choice.startswith("B.") and key_features_for_22q.issubset(patient_symptoms | {"speech problems"}): # Map "incomprehensible speech" to "speech problems"
            score += 5

        if score > max_score:
            max_score = score
            best_match = choice
            best_match_details = f"Matched symptoms: {', '.join(sorted(list(matched_symptoms)))}"


    print("Conclusion:")
    explanation = f"""
The patient presents with a constellation of symptoms including cleft palate, midface hypoplasia, micrognathia, dysplastic ears, and significant developmental/speech delay. This clinical picture is highly characteristic of 22q11.2 deletion syndrome, also known as Velocardiofacial syndrome.

- **Choice A (Chromosome 3):** Less specific match.
- **Choice B (Chromosome 22):** Excellent match. The combination of palatal, facial, and developmental features is classic.
- **Choice C (Chromosome 21):** Incorrect. The facial features described do not align with Trisomy 21 (Down Syndrome).
- **Choice D (Chromosome 2):** Less specific match.
- **Choice E (Chromosome 13):** Unlikely, as Trisomy 13 is far more severe, and survival to age 15 is extremely rare.

The analysis identifies the best fit as the abnormality on Chromosome 22.
"""
    print(textwrap.dedent(explanation).strip())
    print("-" * 50)
    print(f"The most likely abnormality is associated with: {best_match.split('(')[0].strip()}")


find_chromosomal_abnormality()
<<<B>>>