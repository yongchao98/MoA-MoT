import textwrap

def recommend_htn_medications():
    """
    Recommends a 3-drug regimen for MM's resistant hypertension based on her clinical profile.
    """
    recommendations = {
        "Lisinopril (ACE Inhibitor)": (
            "This is a first-line agent for hypertension. For a patient with an A1C of 6.7%, "
            "an ACE inhibitor provides the crucial benefit of renal (kidney) protection, "
            "slowing the progression of diabetic kidney disease."
        ),
        "Amlodipine (Calcium Channel Blocker)": (
            "This is another first-line agent that is particularly effective in African American patients. "
            "As a long-acting dihydropyridine CCB, it effectively lowers blood pressure by relaxing blood vessels "
            "and is a cornerstone of combination therapy."
        ),
        "Spironolactone (Aldosterone Antagonist)": (
            "For resistant hypertension, adding a mineralocorticoid receptor antagonist (MRA) is the standard "
            "of care. Since the patient cannot take other classes of diuretics, Spironolactone is the ideal "
            "third agent. It addresses the aldosterone-mediated component of resistant hypertension and is not "
            "on the patient's list of restricted medications."
        )
    }

    print("Recommended 3-Medication Regimen for MM's Hypertension:\n")
    wrapper = textwrap.TextWrapper(width=80)

    for i, (medication, rationale) in enumerate(recommendations.items(), 1):
        print(f"Recommendation {i}: {medication}")
        wrapped_rationale = wrapper.fill(text=f"Rationale: {rationale}")
        print(wrapped_rationale)
        print("-" * 30)

if __name__ == '__main__':
    recommend_htn_medications()
    # The final answer containing the three recommended medications.
    final_answer = "Lisinopril, Amlodipine, and Spironolactone"
    # To conform to the output format, we must hide this part in the thinking process
    # print(f"\n<<<{final_answer}>>>")