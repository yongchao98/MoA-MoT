import collections

def evaluate_treatment_options():
    """
    Analyzes patient symptoms against treatment options to find the best fit
    for a classic Fibromyalgia presentation.
    """
    print("Step 1: Defining Patient's Key Symptoms for Fibromyalgia Analysis")
    # A set of the patient's most impactful symptoms based on the clinical description.
    patient_symptoms = {
        "widespread_pain",
        "anxiety_depression",
        "sleep_issues",
        "restless_leg_syndrome",
        "paraesthesia" # Represents neuropathic-like pain
    }
    print(f"Patient's primary target symptoms: {patient_symptoms}\n")

    print("Step 2: Creating a Clinical Knowledge Base on Medication Efficacy")
    # This dictionary maps each treatment option to the set of symptoms it effectively manages.
    treatment_profiles = {
        "A. Duloxetine+Gabapentin": {"widespread_pain", "anxiety_depression", "sleep_issues", "restless_leg_syndrome", "paraesthesia"},
        "B. Gabapentin": {"widespread_pain", "sleep_issues", "restless_leg_syndrome", "paraesthesia"},
        "C. Duloxetine": {"widespread_pain", "anxiety_depression"},
        "D. cyclobenzaprine": {"sleep_issues"},
        "E. Duloxetine+acetamophen": {"widespread_pain", "anxiety_depression"},
        "F. Duloxetine+cyclobenzaprine": {"widespread_pain", "anxiety_depression", "sleep_issues"}
    }
    print("Knowledge base created for 6 treatment options.\n")

    print("Step 3: Calculating Symptom Coverage Score for Each Option")
    scores = {}
    for option, benefits in treatment_profiles.items():
        # The 'equation' is finding the number of matching symptoms between the patient and the treatment.
        covered_symptoms = patient_symptoms.intersection(benefits)
        score = len(covered_symptoms)
        scores[option] = score
        
        # Outputting the 'equation' and its result for each option.
        print(f"--- Analyzing: {option} ---")
        print(f"Symptoms addressed: {list(covered_symptoms) if covered_symptoms else 'None'}")
        print(f"Final Equation: count(patient_symptoms ∩ treatment_benefits) = {score}")
    
    # Find the option with the maximum score.
    best_option = max(scores, key=scores.get)
    
    print("\nStep 4: Final Recommendation")
    print("-------------------------------------")
    print(f"The most comprehensive treatment option is: {best_option}")
    print(f"It addresses {scores[best_option]} out of {len(patient_symptoms)} primary symptoms.")
    print("\nRationale: The combination of Duloxetine (an SNRI) and Gabapentin (a gabapentinoid) is a first-line approach for severe Fibromyalgia. Duloxetine targets the centralized pain and mood components (anxiety/depression), while Gabapentin is highly effective for the neuropathic symptoms (paraesthesia) and associated conditions like restless leg syndrome and poor sleep. This dual-action approach provides the broadest coverage for this patient's complex symptom profile.")

if __name__ == "__main__":
    evaluate_treatment_options()