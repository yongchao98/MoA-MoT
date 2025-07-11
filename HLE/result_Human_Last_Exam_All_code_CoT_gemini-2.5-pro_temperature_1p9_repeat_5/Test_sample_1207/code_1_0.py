import textwrap

def solve_clinical_case():
    """
    Analyzes the clinical case and determines the most likely imaging finding.
    """
    
    explanation = """
    The patient presents with a multi-system inflammatory condition affecting the eyes (monocular vision loss), nerves (headaches, hearing loss), joints (pain), lungs (dyspnea), and skin (painful lower extremity lesion). This constellation of symptoms is highly characteristic of Sarcoidosis, a granulomatous disease known for its varied presentations.

    When sarcoidosis affects the central nervous system (Neurosarcoidosis), it can cause cranial nerve palsies (leading to vision and hearing loss) and inflammation of the brain's lining (meninges). The other symptoms also fit well: pulmonary sarcoidosis causes dyspnea, sarcoid arthritis causes joint pain, and erythema nodosum (a common finding in sarcoidosis) causes painful skin nodules on the lower legs.

    Reviewing the answer choices:
    A. Periarticular bone demineralization is for Rheumatoid Arthritis.
    C. Pleural effusion is a non-specific finding.
    D. Vascular hemorrhage is less likely than ischemic events for transient symptoms.
    E. An intrasellar mass would not explain the systemic symptoms.

    B. Leptomeningeal enhancement with "snowball" hyperintensities on MRI is a classic finding for Neurosarcoidosis. This directly corresponds to the most probable underlying diagnosis that unifies all of the patient's symptoms.
    """
    
    final_answer = "B"

    # Use textwrap to format the explanation nicely
    wrapped_explanation = textwrap.fill(textwrap.dedent(explanation).strip(), width=100)
    
    print("Clinical Reasoning:")
    print("-" * 20)
    print(wrapped_explanation)
    print("-" * 20)
    print(f"The most likely image modality and finding is described in choice: {final_answer}")


solve_clinical_case()
print("<<<B>>>")