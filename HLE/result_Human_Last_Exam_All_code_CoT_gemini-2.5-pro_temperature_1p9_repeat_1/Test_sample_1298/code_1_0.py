import textwrap

def solve_medical_case():
    """
    Analyzes the medical case and determines the next best step in management.
    """
    
    print("Step 1: Analyzing the Patient's Presentation")
    analysis_points = [
        "Patient: 20-year-old woman.",
        "Primary Lesion: A mole on her back with classic signs of malignancy (darkening, increasing size, irregular border, raised), highly suspicious for malignant melanoma.",
        "Signs of Metastasis (spread): New dark spots (skin metastases), hip ache (bone metastasis), and cardiopulmonary symptoms (shortness of breath, chest discomfort).",
        "Acute Crisis: The physical exam (muffled heart sounds, jugular venous distention) points to cardiac tamponade, a life-threatening compression of the heart.",
        "Definitive Diagnosis: Pericardiocentesis (draining fluid from around the heart) confirms a malignant pericardial effusion, meaning cancer cells have spread to the pericardium and caused the fluid buildup.",
        "Conclusion from Analysis: The patient has Stage IV (metastatic) malignant melanoma."
    ]
    for point in analysis_points:
        print(textwrap.fill(f"- {point}", width=90))

    print("\nStep 2: Evaluating the Management Options")
    options_evaluation = {
        'A. Prescribe meloxicam': "Incorrect. This is an NSAID for pain/inflammation, it does not treat the underlying widespread cancer.",
        'B. Prescribe low-dose analgesia': "Incorrect. This is for palliative pain control, not the primary treatment for the disease.",
        'D. Chemotherapy to kill the malignant cells': "Correct. The patient has confirmed systemic (metastatic) malignancy. The next crucial step after stabilizing the immediate life threat (tamponade) is to treat the underlying cancer. Chemotherapy is a form of systemic therapy designed to kill cancer cells throughout the body.",
        'E. Immunosuppression': "Incorrect. This would weaken the body's ability to fight cancer. The opposite, immunotherapy, is often used for melanoma.",
        'F. Rapid antibiotic infusion': "Incorrect. There are no signs of a bacterial infection.",
        'G. Radiotherapy to treat the malignant cells': "Incorrect. Radiotherapy is a localized treatment. It is not the primary choice for disease that has spread systemically.",
        'H. Diuretics to reduce the fluid overload': "Incorrect. This is a symptomatic treatment and does not address the cancer causing the fluid. The immediate, life-threatening fluid around the heart has already been addressed by pericardiocentesis."
    }
    
    for option, explanation in options_evaluation.items():
        print(textwrap.fill(f"- {option}: {explanation}", width=90))
        
    print("\nStep 3: Final Conclusion")
    print("The patient's condition is caused by metastatic cancer. Therefore, the next best step is to initiate systemic therapy to treat the cancer throughout her body. Of the choices provided, chemotherapy is the most appropriate option.")

if __name__ == "__main__":
    solve_medical_case()