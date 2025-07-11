import textwrap

def explain_diagnosis():
    """
    Provides a step-by-step explanation for the medical diagnosis based on the clinical vignette.
    """

    patient_profile = {
        "Complaint": "Transient monocular vision loss, pulsatile headaches, joint pain, dyspnea, hearing loss, and a painful lower extremity lesion.",
        "Diagnosis": "The combination of multi-organ symptoms (eye, nerve, lung, joint, skin) is highly suggestive of a systemic granulomatous disease, most notably Sarcoidosis with neurological involvement (Neurosarcoidosis).",
        "Imaging Rationale": "We must identify the imaging finding most characteristic of this diagnosis.",
        "Analysis of Choice B": "'Leptomeningeal enhancement with \"snowball\" hyperintensities visualized by MRI' is a classic and highly specific finding for Neurosarcoidosis.",
        "Final Conclusion": "This MRI finding directly explains the patient's headaches, vision loss, and hearing loss, making it the most expected result given the clinical presentation."
    }

    print("Step-by-Step Reasoning:")
    print("-" * 30)

    for step, explanation in patient_profile.items():
        # Wrap text for better readability in the terminal
        wrapped_text = textwrap.fill(explanation, width=70)
        print(f"Step: {step}\n{wrapped_text}\n")

    final_answer = "B"
    print(f"The correct option is therefore {final_answer}.")

explain_diagnosis()