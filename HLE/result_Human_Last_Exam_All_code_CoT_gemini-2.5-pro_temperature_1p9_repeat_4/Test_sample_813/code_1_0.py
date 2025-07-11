def solve_medical_case():
    """
    This function analyzes the patient's case step-by-step to find the root cause.
    """

    # Patient information
    patient_age = 59
    work_history_duration = 30
    work_industry = "metal smelting, construction, farming"
    family_history = "mood disorders"
    initial_symptoms = ["agitation", "difficulty falling asleep", "increase in sexual activities", "careless spending"]
    later_symptom = "decreased interest in having sex"

    print("Step 1: Analyze the initial symptoms and family history.")
    print(f"The patient's symptoms {initial_symptoms} are characteristic of a manic episode.")
    print(f"The family history of '{family_history}' strongly suggests a predisposition to Bipolar Disorder.")
    print("-" * 30)

    print("Step 2: Determine the likely prescribed medication.")
    print("For a manic episode in Bipolar Disorder, Lithium is a first-line mood-stabilizing medication.")
    print("-" * 30)

    print("Step 3: Analyze the subsequent symptom and connect it to the likely medication.")
    print(f"The patient later developed '{later_symptom}'.")
    print("A primary side effect of Lithium is that it can induce hypothyroidism (underactive thyroid).")
    print("A classic symptom of hypothyroidism is decreased libido or sexual dysfunction.")
    print("-" * 30)

    print("Step 4: Conclude the most probable causal chain.")
    print("The most logical chain of events is:")
    print("  Bipolar Disorder (mania) -> Prescription of Lithium -> Lithium causes Hypothyroidism -> Hypothyroidism causes sexual dysfunction.")
    print("This makes 'A. Lithium induced hypothyroidism' the most comprehensive explanation.")
    print("-" * 30)

# Execute the analysis
solve_medical_case()