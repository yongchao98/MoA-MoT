import sys

def solve_clinical_case():
    """
    Analyzes a clinical case to determine the most likely disease.
    This function will print out the step-by-step reasoning.
    """

    # Step 1: Analyze Patient's Risk Factors
    print("Step 1: Analyzing Patient's Risk Factors")
    age = 62
    smoking_history = "20-pack-year"
    occupation = "ship building"
    print(f" - The patient is a {age}-year-old man.")
    print(f" - He has a significant {smoking_history} smoking history, a major risk factor for lung cancer.")
    print(f" - His work in '{occupation}' suggests potential exposure to asbestos, another major risk factor for lung cancer and mesothelioma.\n")

    # Step 2: Evaluate the Initial Presentation and Underlying Condition
    print("Step 2: Evaluating the Initial Presentation")
    initial_symptoms = ["fatigue", "swelling and pain in wrists, ankles, and elbows (polyarthritis)"]
    progressive_symptoms = ["dizziness", "confusion", "bruising", "difficulty swallowing", "loss of appetite", "shortness of breath"]
    xray_finding = "multiple pulmonary nodules"
    print(f" - Initial symptoms included {', '.join(initial_symptoms)}.")
    print(f" - The condition progressed to include systemic and neurological symptoms: {', '.join(progressive_symptoms)}.")
    print(f" - A chest X-ray revealed '{xray_finding}'.")
    print(" - Conclusion for this step: The combination of risk factors, constitutional symptoms, polyarthritis (a potential paraneoplastic syndrome), and multiple lung nodules strongly suggests an underlying malignancy, most likely lung cancer. The use of steroids for his initial symptoms would have made him immunocompromised.\n")

    # Step 3: Analyze the Acute Illness
    print("Step 3: Analyzing the Acute Illness After Travel")
    acute_symptoms = ["fever", "shortness of breath", "productive cough with green sputum", "cutaneous lesions"]
    treatment = "Aminoglycoside therapy"
    treatment_outcome = "ineffective"
    final_outcome = "died from septic shock"
    print(f" - After traveling, he developed an acute illness with the following symptoms: {', '.join(acute_symptoms)}.")
    print(f" - This presentation is classic for a severe pulmonary infection that has disseminated (spread) to the skin.")
    print(f" - He was treated with '{treatment}', which was '{treatment_outcome}'. Aminoglycosides are not effective against certain types of organisms, such as fungi or weakly acid-fast bacteria like Nocardia.")
    print(f" - His death from '{final_outcome}' indicates the infection was severe, uncontrolled, and systemic.\n")

    # Step 4: Synthesize and Final Diagnosis
    print("Step 4: Synthesizing All Information for a Final Diagnosis")
    print(" - We have an immunocompromised host (due to underlying cancer and steroid use) who develops a severe, disseminated infection.")
    print(" - The infection involves the lungs and skin, and is resistant to aminoglycoside antibiotics.")
    print(" - This clinical picture is highly characteristic of Nocardiosis. Nocardia is an opportunistic bacterium that thrives in immunocompromised individuals. It typically causes pulmonary disease that can spread to the brain (explaining the confusion) and skin (explaining the cutaneous lesions).")
    print(" - The initial underlying disease was likely lung cancer, which created the immunocompromised state necessary for the fatal Nocardia infection to take hold.")

solve_clinical_case()