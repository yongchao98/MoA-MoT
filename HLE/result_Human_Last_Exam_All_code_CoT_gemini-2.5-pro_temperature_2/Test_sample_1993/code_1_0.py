import sys

def find_molecular_abnormality():
    """
    Analyzes a patient's clinical features against a knowledge base of genetic
    disorders to identify the most likely molecular abnormality.
    """
    # Step 1: Define the patient's clinical and genetic features
    patient_features = {
        "ovarian_dysgenesis",
        "short_stature",
        "exercise_intolerance",
        "hypertension",
        "normal_karyotype"
    }

    # Step 2: Create a knowledge base of relevant disorders
    knowledge_base = {
        "Turner Syndrome": {
            "positive_features": {"ovarian_dysgenesis", "short_stature", "hypertension"},
            "negative_features": {"normal_karyotype"},  # Patient's normal karyotype contradicts the typical 45,X
            "abnormality": "Monosomy X (45,X Karyotype)"
        },
        "Mitochondrial Disease": {
            "positive_features": {"ovarian_dysgenesis", "short_stature", "exercise_intolerance", "hypertension", "normal_karyotype"},
            "negative_features": set(),
            "abnormality": "Mutation in the mitochondrial DNA"
        },
        "Swyer Syndrome": {
            "positive_features": {"ovarian_dysgenesis"},
            "negative_features": {"normal_karyotype"}, # Patient is 46,XX; Swyer is 46,XY
            "abnormality": "Mutation in a sex-determining gene on the Y chromosome (e.g., SRY) with a 46,XY karyotype"
        },
        "Isolated Premature Ovarian Insufficiency": {
            "positive_features": {"ovarian_dysgenesis", "normal_karyotype"},
            "negative_features": set(),
            "abnormality": "Mutation in a gene related to ovarian function (e.g., FMR1, BMP15)"
        }
    }

    best_match = None
    max_score = -1

    # Step 3: Implement the matching logic
    for disease, data in knowledge_base.items():
        # Check for contradictions
        if not patient_features.isdisjoint(data["negative_features"]):
            continue  # This diagnosis is disqualified

        # Calculate score based on matching positive features
        score = len(patient_features.intersection(data["positive_features"]))
        
        if score > max_score:
            max_score = score
            best_match = disease

    # Step 4: Print the result
    if best_match:
        print("Based on the clinical presentation, the most likely molecular abnormality is:")
        print(knowledge_base[best_match]["abnormality"])
    else:
        print("Could not determine a likely diagnosis based on the provided features.")

find_molecular_abnormality()