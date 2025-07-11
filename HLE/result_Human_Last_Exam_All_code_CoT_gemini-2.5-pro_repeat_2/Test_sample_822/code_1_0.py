import textwrap

def provide_diagnosis():
    """
    Analyzes the patient's case study and provides the most likely diagnosis.
    """
    reasoning = """
    1. Patient History: The patient is a 62-year-old smoker (20-pack-years) who worked in shipbuilding. This history presents significant risk factors for asbestos exposure and lung cancer.

    2. Initial Symptoms: The initial presentation of pain and swelling in the wrists, ankles, and elbows is characteristic of Hypertrophic Osteoarthropathy (HPOA). HPOA is a well-known paraneoplastic syndrome strongly associated with lung cancer.

    3. Progressive Symptoms & Findings:
       - Chest X-ray: Revealed multiple pulmonary nodules, a common finding in metastatic or advanced lung cancer.
       - Neurological Symptoms: Confusion and dizziness are likely due to another paraneoplastic syndrome, the Syndrome of Inappropriate Antidiuretic Hormone (SIADH), which is particularly common with Small Cell Lung Cancer (SCLC).
       - Other Symptoms: Bruising, difficulty swallowing, and loss of appetite are all consistent with advanced malignancy.

    4. Final Illness: The patient was likely immunocompromised due to the underlying cancer and steroid treatment. The final septic episode was likely an opportunistic infection that did not respond to standard therapy, a common complication in such patients.

    5. Conclusion: The combination of the risk factors, the specific paraneoplastic syndromes (HPOA and SIADH), and the radiological findings of multiple pulmonary nodules makes Small Cell Lung Cancer the most probable underlying disease that explains the entire clinical course.
    """

    final_disease = "Small Cell Lung Cancer with Paraneoplastic Syndromes"

    print("Diagnostic Reasoning:")
    print("--------------------")
    print(textwrap.dedent(reasoning).strip())
    print("\n")
    print("Likely Disease:")
    print("---------------")
    print(final_disease)


provide_diagnosis()

# The final answer is derived from the conclusion of the analysis.
# The patient's entire clinical picture, from the initial joint pain (Hypertrophic Osteoarthropathy)
# to the neurological symptoms (likely SIADH) and lung nodules, points to a single underlying cause.
print("\n<<<Small Cell Lung Cancer with Paraneoplastic Syndromes>>>")