import textwrap

def medical_case_analysis():
    """
    Analyzes a clinical case to determine the most indicative lab parameter
    for a specific condition.
    """
    
    # Case Summary: A 43-year-old female with a history suggestive of Systemic
    # Lupus Erythematosus (SLE) experiences a rapid decline in kidney function
    # after discontinuing corticosteroid treatment.
    
    # Analysis Plan:
    # 1. The patient's history (facial rash, joint pain, hematuria) strongly suggests SLE.
    # 2. The rapid deterioration of kidney function points to a severe flare of Lupus Nephritis.
    # 3. Lupus Nephritis is caused by immune complexes (antibody-antigen clusters) depositing in the kidneys.
    # 4. These immune complexes activate the complement system, a part of the immune response.
    # 5. This activation *consumes* complement proteins, so their levels in the blood drop.
    # 6. Therefore, measuring complement levels provides a direct insight into the cause of the kidney damage.

    print("Step 1: Analyzing the patient's symptoms and history.")
    print("Conclusion: The clinical picture is highly consistent with Systemic Lupus Erythematosus (SLE).")
    print("-" * 50)
    
    print("Step 2: Evaluating the acute event.")
    print("Conclusion: The rapid decline in renal function is characteristic of a severe Lupus Nephritis flare.")
    print("-" * 50)

    print("Step 3: Identifying the key lab indicator for the cause of the flare.")
    
    # Lab parameters and their significance in this context
    # While rising creatinine and anti-dsDNA are important, complement consumption is
    # the most direct marker of the immunologic process causing the damage.
    key_indicator = {
        "Parameter": "Complement Levels",
        "Components": ["C3", "C4"],
        "Expected Finding": "Decreased",
        "Significance": "Indicates active consumption by immune complexes being deposited in the kidneys, which is the direct cause of Lupus Nephritis."
    }
    
    print("The lab parameter that could have best indicated the CAUSE of the rapid renal function decline is:")
    print(f"-> {key_indicator['Parameter']}")

    explanation = (
        f"Specifically, low levels of its components, component 3 ({key_indicator['Components'][0]}) and component 4 ({key_indicator['Components'][1]}), "
        f"would be expected. A sharp decrease in these values signifies that the immune system "
        f"is mounting a significant attack on the kidneys. This process of complement consumption is the direct pathogenic mechanism of a Lupus Nephritis flare."
    )
    
    print("\nExplanation:")
    # Wrap text for better readability in the terminal
    wrapped_explanation = textwrap.fill(explanation, width=70)
    print(wrapped_explanation)

# Execute the analysis
medical_case_analysis()