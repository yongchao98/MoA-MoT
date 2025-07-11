def find_best_lab_indicator():
    """
    Analyzes the clinical scenario to identify the best lab parameter
    indicating the cause of rapid renal decline.
    """
    
    # This dictionary describes relevant lab tests and their significance in this clinical context.
    lab_parameters_info = {
        "Serum Creatinine": "This measures kidney function. It would be elevated, indicating the *result* of the kidney damage, not the primary immunological cause.",
        "Anti-dsDNA Antibodies": "A specific antibody for SLE. High levels indicate high disease activity and can precede a flare, making it a strong predictive marker.",
        "Complement Levels (C3 & C4)": "These immune system proteins are consumed during an active flare, especially in lupus nephritis. Low levels directly reflect the ongoing immune complex deposition and inflammation *causing* the kidney destruction."
    }

    # The patient's history suggests SLE, and the acute symptoms point to a severe lupus nephritis flare.
    # We are looking for the marker that best explains the *cause* of the kidney damage.
    # The consumption of complement is a direct result of the immune process attacking the kidneys.
    best_indicator = "Complement Levels (C3 & C4)"
    explanation = lab_parameters_info[best_indicator]

    print("Analyzing the patient's case to find the best lab indicator for the CAUSE of renal failure...")
    print("-" * 60)
    print(f"The patient's symptoms are highly suggestive of a severe flare of lupus nephritis, which is an autoimmune disease.")
    print(f"The best lab parameter to indicate the underlying immunological cause of this condition is:")
    print(f"\nLab Parameter: {best_indicator}")
    print(f"Reasoning: {explanation}")
    print("-" * 60)
    print("Low C3 and C4 levels would have provided the strongest evidence that an active, complement-consuming autoimmune process was responsible for the rapid decline in the patient's kidney function.")

find_best_lab_indicator()