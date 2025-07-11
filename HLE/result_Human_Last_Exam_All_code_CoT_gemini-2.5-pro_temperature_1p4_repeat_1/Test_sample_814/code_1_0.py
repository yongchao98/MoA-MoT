import collections

def solve_medical_case():
    """
    Analyzes a patient case for Fibromyalgia and determines the best treatment option
    from a given list by scoring each option based on symptom coverage.
    """
    # Step 1: Define the patient's key symptoms that require treatment.
    # The diagnosis is strongly suggested to be Fibromyalgia.
    patient_symptoms = {
        "widespread_pain",
        "anxiety_depression",
        "sleep_issues",
        "restless_leg_syndrome",
        "paresthesia"  # A form of neuropathic pain
    }
    print("Patient Profile Analysis")
    print("="*30)
    print(f"The patient presents with a classic symptom cluster for Fibromyalgia, including:")
    for symptom in sorted(list(patient_symptoms)):
        print(f"- {symptom.replace('_', ' ').title()}")
    print("\n")

    # Step 2: Define the therapeutic targets of each drug for Fibromyalgia.
    drug_efficacy = {
        "Duloxetine": {"widespread_pain", "anxiety_depression", "paresthesia"},
        "Gabapentin": {"widespread_pain", "paresthesia", "sleep_issues", "restless_leg_syndrome"},
        "Cyclobenzaprine": {"sleep_issues"},
        "Acetaminophen": set(),  # Not a primary or effective treatment for fibromyalgia's central pain
        "Ibuprofen": {"widespread_pain"} # Already being taken with minimal effect
    }

    # Step 3: Define the answer choices.
    answer_choices = {
        "A": ["Duloxetine", "Gabapentin"],
        "B": ["Gabapentin"],
        "C": ["Duloxetine"],
        "D": ["Cyclobenzaprine"],
        "E": ["Duloxetine", "Acetaminophen"],
        "F": ["Duloxetine", "Cyclobenzaprine"]
    }

    # Step 4: Evaluate and score each treatment option.
    print("Treatment Option Evaluation")
    print("="*30)
    results = {}
    best_option = ''
    max_score = -1

    for option, drugs in answer_choices.items():
        covered_symptoms = set()
        for drug in drugs:
            covered_symptoms.update(drug_efficacy.get(drug, set()))

        score = len(covered_symptoms)
        results[option] = {
            "score": score,
            "covered": covered_symptoms,
            "drugs": drugs
        }
        if score > max_score:
            max_score = score
            best_option = option

    # Step 5: Print the analysis and conclusion.
    for option, data in sorted(results.items()):
        print(f"Option {option}: {' + '.join(data['drugs'])}")
        print(f"  - Symptom Coverage Score: {data['score']}")
        if data['covered']:
            print(f"  - Covers: {', '.join(sorted(list(data['covered'])))}")
        else:
            print("  - Covers: No primary fibromyalgia symptoms")
        print("-" * 25)

    print("\nConclusion")
    print("="*30)
    print("The ideal treatment should be comprehensive.")
    print(f"Option {best_option} achieves the highest symptom coverage score of {max_score}.")
    print("It combines Duloxetine to treat pain and depression with Gabapentin, which targets a different pain pathway and also specifically addresses sleep issues and restless leg syndrome.")
    print("This combination provides the most robust and synergistic therapy for this patient's full range of symptoms.")
    print(f"\nThe best choice is: {best_option}")

solve_medical_case()