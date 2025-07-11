import sys

def solve_clinical_case():
    """
    This function analyzes a clinical vignette to determine the expected location of a rash.
    It prints the step-by-step reasoning and the final conclusion.
    """
    print("Step 1: Analyze the key clinical findings.")
    key_symptoms = "muscle weakness, myalgia (muscle pain), arthralgia (joint pain), and fatigue."
    key_physical_finding = "periorbital erythema (redness around the eyes)."
    print(f"The patient presents with systemic symptoms like {key_symptoms}")
    print(f"The key physical finding is {key_physical_finding}")
    print("-" * 30)

    print("Step 2: Synthesize the findings into a likely diagnosis.")
    print("The combination of muscle inflammation (myositis) and skin manifestations (dermatitis) points strongly towards Dermatomyositis.")
    print("-" * 30)

    print("Step 3: Relate the diagnosis to the specific physical finding.")
    print("In Dermatomyositis, 'periorbital erythema' is the clinical description for a pathognomonic (uniquely characteristic) rash called the 'heliotrope rash'.")
    print("A heliotrope rash is a violaceous or reddish rash that appears on and around the eyelids.")
    print("-" * 30)

    print("Step 4: Determine the anatomical location.")
    print("The physical exam finding of 'periorbital erythema' directly describes a rash located on the anatomical region of the eyelids.")
    print("Therefore, a rash is not just expected on the eyelids, it is already described there in the vignette.")
    print("-" * 30)
    
    answer_choice = "C"
    answer_text = "Eyelids"
    print(f"Conclusion: The correct answer is Choice {answer_choice}, which corresponds to the anatomical region: {answer_text}.")

# Execute the analysis
solve_clinical_case()