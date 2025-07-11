def recommend_htn_medications():
    """
    Recommends a three-drug regimen for the patient MM based on her clinical profile.
    """
    # Patient Profile Key Points:
    # - Resistant Stage II Hypertension (BP 145-153/85-89)
    # - African American
    # - Comorbidities: Type 2 Diabetes (A1C 6.7%), Tachycardia (HR 91)
    # - Restrictions: Unable to take common diuretics, MRAs, and verapamil.

    # Recommended Drug 1: Angiotensin Receptor Blocker (ARB)
    drug_1 = {
        "class": "Angiotensin Receptor Blocker (ARB)",
        "example": "Losartan",
        "rationale": ("Provides guideline-directed treatment for hypertension in a patient with diabetes, "
                      "offering crucial renal protection. It works by blocking the renin-angiotensin system.")
    }

    # Recommended Drug 2: Dihydropyridine Calcium Channel Blocker (CCB)
    drug_2 = {
        "class": "Dihydropyridine Calcium Channel Blocker (CCB)",
        "example": "Amlodipine",
        "rationale": ("A first-line agent that is particularly effective in African American patients. "
                      "It lowers blood pressure through vasodilation, providing a synergistic effect with the ARB.")
    }

    # Recommended Drug 3: Beta-1 Selective Blocker
    drug_3 = {
        "class": "Beta-1 Selective Blocker",
        "example": "Metoprolol Succinate",
        "rationale": ("Addresses the patient's elevated heart rate (tachycardia), which contributes to her high blood pressure. "
                      "This class is a logical third agent when diuretics and MRAs are contraindicated.")
    }

    print("Recommendation for MM's Hypertension Treatment:")
    print("-" * 50)
    
    print(f"1. Medication Class: {drug_1['class']}")
    print(f"   Recommended Drug: {drug_1['example']}")
    print(f"   Rationale: {drug_1['rationale']}\n")

    print(f"2. Medication Class: {drug_2['class']}")
    print(f"   Recommended Drug: {drug_2['example']}")
    print(f"   Rationale: {drug_2['rationale']}\n")

    print(f"3. Medication Class: {drug_3['class']}")
    print(f"   Recommended Drug: {drug_3['example']}")
    print(f"   Rationale: {drug_3['rationale']}\n")

    print("-" * 50)
    print("This three-drug regimen provides a multi-mechanistic approach to control her resistant hypertension while addressing her specific comorbidities and adhering to her medication restrictions.")

if __name__ == '__main__':
    recommend_htn_medications()