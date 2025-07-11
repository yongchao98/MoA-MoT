import textwrap

def analyze_patient_case():
    """
    Analyzes a clinical case to determine the significance of a new diet.
    """

    # --- Patient & Medical Information ---
    patient_condition = "Schizophrenia"
    medication_mechanism = "Dopamine D2 receptor antagonism"
    
    # --- Diet Information ---
    # "Bean salad" is a clue pointing to beans rich in L-DOPA.
    food_description = "Diet that tastes like bean salad"
    likely_food_source = "Fava beans or Velvet beans"
    active_compound_in_food = "Levodopa (L-DOPA)"
    compound_effect = "Metabolic precursor to Dopamine"

    # --- Analysis Logic ---
    print("Patient Medical Analysis:")
    print("-" * 30)
    print(f"Initial Diagnosis: {patient_condition}")
    print(f"Probable Postpartum Diagnosis: Sheehan's Syndrome (based on symptoms)")
    print(f"Primary Medication Action: {medication_mechanism}")
    print("\nDietary Analysis:")
    print("-" * 30)
    print(f"Diet Clue: {food_description}")
    print(f"Likely Ingredient: {likely_food_source}")
    print(f"Key Substance in Ingredient: {active_compound_in_food}")
    print(f"Biochemical Effect: This substance is a {compound_effect}.\n")

    # --- Conclusion ---
    print("Clinical Conclusion and Importance of the Food:")
    print("-" * 50)
    
    explanation = (
        f"The patient's {patient_condition} is treated with medication that blocks dopamine. "
        f"The new food, likely containing sources like {likely_food_source}, is important because "
        f"it is a rich source of {active_compound_in_food}. L-DOPA is converted into dopamine in the body. "
        "Consuming this food introduces a high amount of a dopamine precursor, which directly "
        "counteracts the effect of her antipsychotic medication. This creates a significant risk of "
        "triggering a severe relapse of her schizophrenic symptoms, such as false beliefs and visions."
    )
    
    # Use textwrap for clean output
    wrapped_explanation = textwrap.fill(explanation, width=70)
    print(wrapped_explanation)


# Execute the analysis
if __name__ == "__main__":
    analyze_patient_case()