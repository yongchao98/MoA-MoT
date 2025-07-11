def analyze_patient_case():
    """
    Analyzes the patient's symptoms and evaluates the provided anatomical choices.
    """
    symptoms = {
        "Left-sided facial weakness": "Points to a lesion of the Left Facial Nerve (Cranial Nerve VII).",
        "Inability to lift left eyebrow": "Specifically involves the frontalis muscle, innervated by the Facial Nerve (CN VII).",
        "Loss of left acoustic reflex": "The efferent limb of this reflex is controlled by the stapedius muscle, innervated by the Facial Nerve (CN VII).",
        "Hoarseness and cough": "These symptoms point to dysfunction of the larynx, which is primarily innervated by the Vagus Nerve (Cranial Nerve X), especially its recurrent laryngeal branch.",
        "Small mass in thoracic cavity": "This is a critical finding. The left recurrent laryngeal nerve (a branch of the Vagus Nerve) takes a long path, looping under the aorta in the thorax. A mass here can easily compress this nerve, causing vocal cord paralysis and hoarseness."
    }

    print("Step 1: Analyzing the Patient's Symptoms and Findings\n" + "="*55)
    for symptom, analysis in symptoms.items():
        print(f"- {symptom}: {analysis}")

    print("\nStep 2: Synthesizing the Information\n" + "="*55)
    print("The patient presents with signs of both Facial Nerve (CN VII) and Vagus Nerve (CN X) palsy on the left side.")
    print("The hoarseness linked with a thoracic mass is a classic presentation for compression of the recurrent laryngeal nerve.")
    print("This nerve controls most of the muscles for voice production in the larynx.\n")

    print("Step 3: Evaluating the Anatomical Answer Choices\n" + "="*55)
    choices = {
        "A. Tensor tympani": "Innervated by CN V (Trigeminal). While involved in hearing, the primary evidence points to CN VII and X.",
        "B. Lateral rectus": "Innervated by CN VI (Abducens). The patient has no reported eye movement issues.",
        "C. Intercostal muscles": "Innervated by intercostal nerves. Not directly related to the primary cranial nerve symptoms.",
        "D. Cricothyroid": "A key muscle of the larynx, innervated by a branch of the Vagus Nerve (CN X). Its function is critical for voice pitch. Dysfunction of the laryngeal muscles, like the cricothyroid, directly causes hoarseness. This choice directly relates a key symptom (hoarseness) to the nerve (Vagus) implicated by the thoracic mass.",
        "E. Stylopharyngeus": "Innervated by CN IX (Glossopharyngeal). The symptoms do not align with a primary CN IX lesion."
    }

    for choice, explanation in choices.items():
        print(f"- {choice}: {explanation}")

    print("\nStep 4: Conclusion\n" + "="*55)
    print("The cricothyroid muscle is part of the larynx, the organ responsible for voice.")
    print("The symptom of hoarseness is a direct result of laryngeal muscle paralysis.")
    print("This paralysis is best explained by the thoracic mass compressing the vagus nerve or its recurrent laryngeal branch.")
    print("Therefore, the cricothyroid is the most important anatomical structure to consider among the choices, as it directly connects a major symptom (hoarseness) to the underlying pathology (thoracic mass affecting the vagus nerve).")


if __name__ == "__main__":
    analyze_patient_case()