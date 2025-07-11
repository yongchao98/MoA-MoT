def find_most_likely_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by evaluating each option against the key findings.
    """

    # --- Step 1: Define the key clinical findings from the case ---
    print("Step 1: Summarizing the key clinical findings.")
    findings = {
        "Procedure": "Difficult colonoscopy with no polypectomy.",
        "Pain Location": "Left Upper Quadrant (LUQ) and left shoulder (Kehr's sign).",
        "Hemodynamics": "Acute hypovolemic shock (tachycardia, hypotension).",
        "Key Lab Change": "Massive drop in hemoglobin, indicating severe hemorrhage.",
        "Physical Exam": "Abdominal distension and peritoneal signs."
    }
    for key, value in findings.items():
        print(f"- {key}: {value}")

    print("\nStep 2: Evaluating each potential diagnosis.")

    # --- Step 2: Analyze each diagnosis ---

    # A. Colonic Perforation
    print("\nAnalysis of A. Colonic Perforation:")
    print(" - Consistent with post-colonoscopy complication: Yes.")
    print(" - Explains peritoneal signs: Yes.")
    print(" - Explains severe hemorrhage: Possible, but less common.")
    print(" - Explains LUQ + shoulder pain specifically: No. This pain pattern is not typical for colonic perforation.")
    print("   -> Verdict: Poor fit for the specific pain location.")

    # B. Lower GI Bleeding
    print("\nAnalysis of B. Lower GI Bleeding:")
    print(" - Consistent with post-colonoscopy complication: Yes.")
    print(" - Explains severe hemorrhage: Yes.")
    print(" - Explains LUQ + shoulder pain: No. Pain is typically crampy/lower abdominal, not LUQ with shoulder radiation.")
    print("   -> Verdict: Inconsistent with the location of pain.")

    # C. Splenic Laceration
    print("\nAnalysis of C. Splenic Laceration:")
    print(" - Consistent with difficult colonoscopy (traction injury): Yes.")
    print(" - Explains LUQ pain (location of spleen): Yes, perfect match.")
    print(" - Explains left shoulder pain (Kehr's sign from diaphragmatic irritation): Yes, classic sign.")
    print(" - Explains severe hemorrhage and shock (spleen is highly vascular): Yes, perfect match.")
    print("   -> Verdict: Excellent fit. Explains all key findings.")

    # D. Postpolypectomy Syndrome
    print("\nAnalysis of D. Postpolypectomy Syndrome:")
    print(" - Requires a polypectomy to have been performed.")
    print(" - Case states: 'No polypectomy was performed'.")
    print("   -> Verdict: Ruled out by the case history.")

    # --- Step 3: Present the final conclusion with numerical evidence ---
    print("\nStep 3: Concluding based on the analysis.")
    print("The diagnosis of Splenic Laceration is the only one that cohesively explains the entire clinical picture.")
    print("Let's review the key numbers that point to massive hemorrhage:")

    initial_hgb = 11.7
    later_hgb = 6.5
    hgb_drop = initial_hgb - later_hgb

    print(f"The patient's hemoglobin level dropped from {initial_hgb} g/dL to {later_hgb} g/dL.")
    print("Final Equation for Hemoglobin Drop:")
    print(f"{initial_hgb} - {later_hgb} = {hgb_drop:.1f} g/dL")
    print("This significant drop, combined with the classic pain pattern, strongly supports a diagnosis of splenic laceration and associated hemoperitoneum.")

    # --- Final Answer ---
    final_answer = "C"
    print(f"\nThe most likely diagnosis is C.")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    find_most_likely_diagnosis()