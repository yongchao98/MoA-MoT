import pandas as pd

def find_best_treatment():
    """
    Analyzes a clinical case of suspected fibromyalgia and determines the best treatment
    option by scoring each choice against the patient's main symptoms.
    """
    # Step 1: Define patient's key symptoms and drug efficacy profiles
    symptoms = ["Widespread Pain", "Anxiety/Depression", "Sleep Issues", "Neuropathic Symptoms"]
    
    # Efficacy scores (0=none, 1=low, 2=moderate, 3=high)
    drug_efficacy = {
        "Duloxetine":         {"Widespread Pain": 3, "Anxiety/Depression": 3, "Sleep Issues": 1, "Neuropathic Symptoms": 1},
        "Gabapentin":         {"Widespread Pain": 2, "Anxiety/Depression": 0, "Sleep Issues": 2, "Neuropathic Symptoms": 3},
        "Cyclobenzaprine":    {"Widespread Pain": 1, "Anxiety/Depression": 0, "Sleep Issues": 2, "Neuropathic Symptoms": 0},
        "Acetaminophen":      {"Widespread Pain": 0, "Anxiety/Depression": 0, "Sleep Issues": 0, "Neuropathic Symptoms": 0}
    }

    # Step 2: Define the answer choices
    treatment_options = {
        "A": ["Duloxetine", "Gabapentin"],
        "B": ["Gabapentin"],
        "C": ["Duloxetine"],
        "D": ["Cyclobenzaprine"],
        "E": ["Duloxetine", "Acetaminophen"],
        "F": ["Duloxetine", "Cyclobenzaprine"]
    }

    # Step 3 & 4: Calculate scores and find the best option
    results = {}
    print("Evaluating treatment options based on efficacy scores:\n")

    for option, drugs in treatment_options.items():
        total_score = 0
        equation_parts = []
        for drug in drugs:
            drug_score = sum(drug_efficacy[drug].values())
            total_score += drug_score
            equation_parts.append(f"{drug}({drug_score})")
        
        equation = " + ".join(equation_parts)
        results[option] = total_score
        print(f"Option {option}: { ' + '.join(drugs) }")
        print(f"Score Calculation: {equation} = {total_score}\n")

    best_option = max(results, key=results.get)
    
    print("--- Conclusion ---")
    print(f"The patient's symptoms include widespread pain, anxiety/depression, sleep issues, and neuropathic symptoms.")
    print(f"The combination in Option {best_option} provides the most comprehensive coverage for this constellation of symptoms.")
    print(f"Therefore, {best_option} is the best choice.")

# Run the analysis
find_best_treatment()