import textwrap

def explain_diagnosis():
    """
    Analyzes the patient's symptoms and explains the most likely diagnosis and corresponding imaging finding.
    """

    # 1. Deconstruct the clinical presentation
    age = 44
    symptoms = {
        "Ophthalmic": "Transient monocular vision loss",
        "Neurologic": "Pulsatile headaches",
        "Auditory": "Hearing loss",
        "Musculoskeletal": "Joint pain (arthralgia)",
        "Pulmonary": "Dyspnea (shortness of breath)",
        "Dermatologic": "Painful area on lower extremity (suggestive of erythema nodosum)"
    }

    print(f"Patient Profile: A {age}-year-old woman presents with a complex, multi-system disorder.")
    print("-" * 30)

    # 2. Rationale for Diagnosis
    print("Clinical Reasoning:")
    reasoning_text = """
    The combination of symptoms affecting multiple, seemingly unrelated organ systems (eyes, brain, ears, lungs, joints, skin) is highly suggestive of a systemic inflammatory condition. Sarcoidosis is a multisystem granulomatous disease that commonly presents this way.

    - The neurologic, vision, and hearing symptoms (headaches, vision loss, hearing loss) point to Central Nervous System (CNS) involvement, specifically neurosarcoidosis.
    - Neurosarcoidosis often involves inflammation of the cranial nerves (like the optic and vestibulocochlear nerves) and the meninges (the membranes covering the brain).
    - The systemic symptoms like joint pain (sarcoid arthritis), dyspnea (pulmonary sarcoidosis), and painful skin nodules (erythema nodosum) are also classic features of sarcoidosis.
    """
    print(textwrap.fill(reasoning_text, 80))
    print("-" * 30)

    # 3. Evaluate the Imaging Finding
    print("Evaluation of Imaging Choice B:")
    imaging_finding_explanation = """
    The question asks for the expected image modality and finding. Choice B is 'Leptomeningeal enhancement with "snowball" hyperintensities visualized by MRI'.

    - Modality: MRI is the gold standard for visualizing soft tissues of the brain and meninges.
    - Finding: Leptomeningeal enhancement on a contrast-enhanced MRI indicates inflammation of the meninges, which directly explains the patient's headaches and cranial nerve-related symptoms (vision and hearing loss). Parenchymal "snowball" or "cotton ball" lesions are also known findings in neurosarcoidosis. This choice provides a unifying diagnosis for the most specific and concerning symptoms.
    """
    print(textwrap.fill(imaging_finding_explanation, 80))
    print("-" * 30)

    # 4. Final Answer
    final_answer = "B"
    print(f"The most fitting answer is B, as it describes the characteristic MRI findings of neurosarcoidosis, a condition that explains the patient's entire constellation of symptoms.")

    # The final answer in the required format
    print(f"\n<<<B>>>")

explain_diagnosis()