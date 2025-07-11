def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis from a set of choices.
    """
    # Step 1: Define the clinical information
    patient_info = {
        "age": "1-year-old",
        "findings": ["hypertrophic scarring", "erythema", "spasticity"],
        "labs": "Negative for anti-Mi-2"
    }

    print("Analyzing the Clinical Case:")
    print(f"Patient Age: {patient_info['age']}")
    print(f"Key Findings: {', '.join(patient_info['findings'])}")
    print(f"Lab Results: {patient_info['labs']}\n")

    # Step 2: Evaluate each diagnosis
    print("Evaluating the Answer Choices:\n")

    # Choice A: Ectropion
    print("A. Ectropion:")
    print("   - This is an ophthalmologic condition where the eyelid turns outward.")
    print("   - It does not explain systemic findings like hypertrophic scarring on the body, erythema (unless localized to the eye), or spasticity.")
    print("   - Conclusion: Unlikely.\n")

    # Choice B: McArdle disease
    print("B. McArdle disease (Glycogen Storage Disease Type V):")
    print("   - This is a metabolic disorder affecting muscle energy.")
    print("   - Typical onset is in adolescence or adulthood, not infancy.")
    print("   - Symptoms include exercise intolerance, muscle cramps, and fatigue, not spasticity or hypertrophic scarring.")
    print("   - Conclusion: Unlikely.\n")

    # Choice C: Dermatomyositis
    print("C. Dermatomyositis:")
    print("   - This is an inflammatory disease affecting the muscles and skin. The juvenile form (JDM) occurs in children.")
    print("   - Erythema (redness) is a classic skin manifestation (e.g., heliotrope rash, Gottron's papules).")
    print("   - While the primary muscle symptom is weakness, severe, chronic inflammation can lead to joint contractures and muscle shortening, causing stiffness that might be described as 'spasticity'.")
    print("   - Calcinosis cutis (calcium deposits in the skin), a common complication of JDM, can ulcerate and lead to significant scarring, which could be hypertrophic.")
    print("   - The anti-Mi-2 antibody is highly specific but not sensitive; it is negative in a large majority of JDM cases, so a negative result does not rule out the diagnosis.")
    print("   - Conclusion: This is the most plausible diagnosis, accounting for the combination of skin and severe muscle-related findings in a child, despite the atypical terminology.\n")

    # Choice D: McCune Albright syndrome
    print("D. McCune Albright syndrome:")
    print("   - This is a genetic disorder with a classic triad of fibrous dysplasia (bone disease), café-au-lait skin spots, and endocrine problems.")
    print("   - The patient's symptoms of erythema, spasticity, and hypertrophic scarring do not match this triad.")
    print("   - Conclusion: Unlikely.\n")

    # Choice E: Cataracts
    print("E. Cataracts:")
    print("   - This is a clouding of the lens in the eye.")
    print("   - It is an isolated eye condition and does not explain the other systemic symptoms.")
    print("   - Conclusion: Unlikely.\n")
    
    # Step 3: Final Conclusion
    print("-----------------------------------")
    print("Final Determination:")
    print("Dermatomyositis is the best fit. It explains the skin finding (erythema), can explain the hypertrophic scarring as a complication of calcinosis, and the 'spasticity' as a manifestation of severe muscle contractures. The negative anti-Mi-2 lab is expected in many cases of juvenile dermatomyositis.")
    print("The most likely diagnosis is C.")
    print("-----------------------------------")


if __name__ == "__main__":
    solve_medical_case()