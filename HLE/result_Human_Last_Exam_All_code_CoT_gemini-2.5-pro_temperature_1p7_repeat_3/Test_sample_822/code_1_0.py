import sys

def analyze_case():
    """
    Analyzes the patient's case to determine the most likely diagnosis.
    """

    # Patient Data
    age = 62
    pack_years = 20
    occupation = "Ship Building"
    
    print("Patient Profile and Risk Factor Analysis:")
    print(f" - Age: {age} years old.")
    print(f" - Smoking History: {pack_years} pack-years. This is a major risk factor for lung cancer.")
    print(f" - Occupational History: Work in '{occupation}' implies significant exposure to asbestos, which dramatically increases the risk of lung cancer, especially in smokers.")
    print("-" * 40)
    
    print("Clinical Presentation Analysis:")
    print(" - Initial Symptoms: Symmetrical pain and swelling in wrists, ankles, and elbows are characteristic of Hypertrophic Pulmonary Osteoarthropathy (HPOA), a paraneoplastic syndrome strongly associated with non-small cell lung cancer.")
    print(" - Imaging: Chest X-ray showing 'multiple pulmonary nodules' is a classic sign of lung cancer, either as multiple primary tumors or metastases.")
    print(" - Disease Progression: Symptoms like difficulty swallowing, confusion, and shortness of breath indicate advanced, likely metastatic, disease.")
    print(" - Terminal Illness: The development of a severe infection unresponsive to initial antibiotics, leading to septic shock, is a common terminal event in patients who are immunocompromised due to an underlying malignancy and treatments like steroids.")
    print("-" * 40)

    print("Equation of Evidence:")
    print("The combination of key factors leads to the diagnosis.")
    print(f"({age}-year-old man) + ({pack_years} pack-year smoking history) + (Asbestos Exposure) + (Paraneoplastic Syndrome) + (Pulmonary Nodules) => Most Likely Diagnosis")
    print("-" * 40)
    
    final_diagnosis = "Lung Cancer"
    print(f"Conclusion: The combination of profound risk factors (smoking, asbestos) with classic signs (paraneoplastic syndrome of the joints, pulmonary nodules) makes {final_diagnosis} the most probable underlying disease.")

if __name__ == "__main__":
    analyze_case()
    # The final answer is wrapped in the specified format below.
    # We use sys.stdout.write to avoid adding an extra newline before the final delimiter.
    sys.stdout.write("<<<Lung Cancer>>>")
