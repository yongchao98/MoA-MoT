import textwrap

def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Wrap text for cleaner printing
    wrapper = textwrap.TextWrapper(width=80)

    # Patient Profile
    patient_profile = {
        "Condition": "Resistant Hypertension (Stage II)",
        "Goal BP": "< 130/90 mmHg",
        "Key Comorbidities": "Type 2 Diabetes (A1C=6.7%), Hypercholesterolemia, Hypothyroidism",
        "Demographic": "27-year-old African American female"
    }

    # Rationale
    print("Step 1: Rationale for Medication Selection")
    print("------------------------------------------")
    intro_text = (
        "The patient has resistant hypertension, defined as uncontrolled blood pressure despite being on "
        "three medications. She also has Type 2 Diabetes, a compelling indication for specific "
        "antihypertensive classes. The goal is to recommend an optimized, potent, guideline-directed "
        "3-drug regimen that respects her medication intolerance list."
    )
    print("\n".join(wrapper.wrap(intro_text)))
    print("-" * 42)

    # Medication Recommendations
    medications = {
        "Lisinopril": (
            "An ACE (Angiotensin-Converting Enzyme) Inhibitor. This class is strongly recommended as "
            "a first-line agent for patients with hypertension and diabetes, as it provides "
            "significant kidney protection."
        ),
        "Amlodipine": (
            "A Dihydropyridine Calcium Channel Blocker (CCB). CCBs are highly effective in treating "
            "hypertension in African American patients and provide potent blood pressure reduction. "
            "It is not on the patient's exclusion list."
        ),
        "Chlorthalidone": (
            "A thiazide-like diuretic. An effective 3-drug regimen for resistant hypertension typically "
            "includes a diuretic. While several common diuretics are excluded for this patient, "
            "Chlorthalidone is not. It is more potent and longer-acting than hydrochlorothiazide and is "
            "often preferred in clinical guidelines."
        )
    }

    print("\nStep 2: Recommended 3-Drug Regimen")
    print("-----------------------------------")
    for med, reason in medications.items():
        print(f"\nRecommendation 1: {med}")
        print("   Reason: " + "\n".join(wrapper.wrap(reason, initial_indent='   ', subsequent_indent='   ')))
        
    print("\n\nSummary:")
    print("This combination of an ACE inhibitor, a calcium channel blocker, and a potent diuretic creates a synergistic regimen that addresses the patient's hypertension from multiple angles, is consistent with evidence-based guidelines for her specific conditions, and avoids her documented medication intolerances.")

# Execute the function to print the recommendations
recommend_htn_medications()
<<<Lisinopril, Amlodipine, Chlorthalidone>>>