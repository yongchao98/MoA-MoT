import textwrap

def solve_medical_case():
    """
    This function analyzes the clinical vignette step-by-step to arrive at the most logical answer.
    """

    print("Step 1: Analyze the patient's key symptoms and findings.")
    symptoms = {
        "Left-sided facial weakness (cannot lift eyebrow)": "Points to a peripheral lesion of the Facial Nerve (CN VII).",
        "Loss of acoustic reflex in the left ear": "The stapedius muscle, responsible for the acoustic reflex, is innervated by the Facial Nerve (CN VII). This reinforces CN VII involvement.",
        "Hoarseness and cough": "These symptoms strongly suggest an issue with the Vagus Nerve (CN X), specifically its recurrent laryngeal branch which innervates the vocal cords.",
        "Small mass in the thoracic cavity": "This is a critical finding that can explain the Vagus Nerve (CN X) symptoms, as the nerve travels through the thorax."
    }
    for symptom, explanation in symptoms.items():
        print(f"- {symptom}: {explanation}")

    print("\nStep 2: Synthesize the findings.")
    print(textwrap.fill(
        "The patient presents with signs of both CN VII and CN X palsy on the left side. "
        "The hoarseness is a classic sign of recurrent laryngeal nerve palsy. The left recurrent laryngeal nerve has a long course, looping under the aortic arch in the thorax before ascending to the larynx. "
        "A thoracic mass can easily compress this nerve, leading to vocal cord paralysis and hoarseness."
    , 80))

    print("\nStep 3: Evaluate the answer choices based on the synthesis.")
    choices = {
        "A. Tensor tympani": "Innervated by CN V (Trigeminal Nerve). Not consistent with the primary findings.",
        "B. Lateral rectus": "Innervated by CN VI (Abducens Nerve). The patient has no noted eye movement issues.",
        "C. Intercostal muscles": "Involved in breathing, but the key localizing signs are the cranial nerve deficits. These muscles are not the primary structure of concern.",
        "D. Cricothyroid": "This is an intrinsic muscle of the larynx. The patient's hoarseness is a direct result of laryngeal muscle dysfunction. While most laryngeal muscles are innervated by the recurrent laryngeal nerve, the cricothyroid is innervated by the superior laryngeal nerve (also a branch of the Vagus nerve). Because hoarseness is a cardinal symptom pointing to a problem with the larynx, the laryngeal muscles themselves are a central part of the diagnosis. The Cricothyroid is the only laryngeal muscle listed as an option.",
        "E. Stylopharyngeus": "Innervated by CN IX (Glossopharyngeal Nerve). While involved in swallowing, the key symptoms point more directly to CN VII and CN X."
    }
    for choice, explanation in choices.items():
        print(f"- {choice}: {explanation}")

    print("\nStep 4: Conclude the most important structure.")
    print(textwrap.fill(
        "The most prominent and localizing symptom explained by the thoracic mass is hoarseness, which is a disorder of the larynx. Therefore, an anatomical structure within the larynx is the most relevant to consider. Of the choices provided, the Cricothyroid is the only intrinsic muscle of the larynx."
    , 80))

if __name__ == "__main__":
    solve_medical_case()
    print("\nFinal Answer:")
    print("<<<D>>>")
