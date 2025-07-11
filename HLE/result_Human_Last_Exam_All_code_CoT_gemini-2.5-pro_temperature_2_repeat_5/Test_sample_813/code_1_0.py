import textwrap

def analyze_patient_case():
    """
    Analyzes the patient case to determine the root cause of sexual dysfunction.
    """
    # 1. Deconstruct the patient's profile and symptoms.
    initial_symptoms = [
        "agitation",
        "difficulty falling asleep",
        "increase in sexual activities",
        "careless spending"
    ]
    family_history = "mood disorders"
    occupational_history = "30-year work history of metal smelting"
    new_symptom = "decreased interest in having sex (sexual dysfunction)"

    # 2. Print the reasoning process step-by-step.
    print("Step 1: Initial Diagnosis Inference")
    print(textwrap.fill(
        f"The patient's initial symptoms ({', '.join(initial_symptoms)}) are classic indicators of a manic episode. This, combined with a family history of '{family_history}', strongly suggests a diagnosis of Bipolar Disorder.",
        width=80
    ))
    print("-" * 20)

    print("Step 2: Treatment Inference")
    print(textwrap.fill(
        "A standard and very common medication prescribed to treat behavioral disturbances associated with Bipolar Disorder mania is Lithium.",
        width=80
    ))
    print("-" * 20)

    print("Step 3: Connecting Treatment to New Symptom")
    print(textwrap.fill(
        f"The patient later developed a new symptom: '{new_symptom}'. We must evaluate the potential causes based on the likely treatment.",
        width=80
    ))
    print("-" * 20)

    print("Step 4: Evaluating the Most Likely Cause")
    print(textwrap.fill(
        "Lithium is well-known to cause hypothyroidism as a side effect. A primary symptom of hypothyroidism is decreased libido and sexual dysfunction. This creates a direct and logical causal chain: Bipolar Mania -> Lithium Treatment -> Hypothyroidism -> Sexual Dysfunction.",
        width=80
    ))
    print("-" * 20)

    print("Step 5: Conclusion")
    print(textwrap.fill(
        f"While the '{occupational_history}' is a risk factor for heavy metal toxicity (Options B, C, D, E), it does not explain the specific sequence of events as well as Lithium-induced hypothyroidism. The sexual dysfunction appeared *after* the medication was started, pointing to an iatrogenic cause.",
        width=80
    ))
    print("-" * 20)
    print("Final Answer based on the analysis is A.")


if __name__ == '__main__':
    analyze_patient_case()
    print("\n<<<A>>>")
