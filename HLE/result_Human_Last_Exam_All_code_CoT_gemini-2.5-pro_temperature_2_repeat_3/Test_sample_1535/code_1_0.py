import sys

def solve_medical_case():
    """
    This function analyzes the clinical vignette to identify the correct anatomical region for a rash.
    """
    # Patient's key symptoms and signs:
    symptoms = {
        "History": ["myalgia", "muscle weakness", "fatigue", "arthralgia"],
        "Physical Exam": ["periorbital erythema"]
    }

    # Reasoning process:
    # 1. The combination of muscle inflammation (myalgia, weakness) and skin findings points to Dermatomyositis.
    diagnosis = "Dermatomyositis"
    
    # 2. A key pathognomonic (highly specific) skin sign of Dermatomyositis is the Heliotrope rash.
    # 3. The Heliotrope rash is located on the eyelids.
    location_of_heliotrope_rash = "Eyelids"

    # 4. The patient's physical exam finding, "periorbital erythema" (redness around the eyes), is a clinical description of the Heliotrope rash.
    # 5. Therefore, the eyelids are the anatomical region expected to have the characteristic rash.
    
    print("Rationale:")
    print(f"The patient's combination of muscle weakness and skin findings strongly suggests a diagnosis of {diagnosis}.")
    print("A classic and highly specific skin manifestation of this condition is the Heliotrope rash.")
    print(f"This rash is characteristically found on the '{location_of_heliotrope_rash}'.")
    print("The patient's exam notes 'periorbital erythema', which directly corresponds to a rash in this location.")
    
    # Mapping answer choices to potential signs
    answer_choices = {
        'A': "Dorsum of the hands (Gottron's papules)",
        'B': "Nose (part of a malar rash, less specific)",
        'C': "Eyelids (Heliotrope rash)",
        'D': "Groin (not a typical location)",
        'E': "Shoulders (Shawl sign)"
    }
    
    print("\nBased on the analysis, the most specific and expected location for the rash described is C.")
    # The final answer is determined to be 'C' based on the clinical evidence.
    final_answer = 'C'


solve_medical_case()