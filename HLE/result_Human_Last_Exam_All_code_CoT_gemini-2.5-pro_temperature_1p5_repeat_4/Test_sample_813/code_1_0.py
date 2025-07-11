import textwrap

def solve_clinical_case():
    """
    This function analyzes the clinical vignette step-by-step to determine the root cause.
    """

    # Step 1: Analyze the patient's initial symptoms.
    print("Step 1: Analyzing the patient's initial presentation.")
    symptoms = "Agitation, difficulty falling asleep, an increase in sexual activities, and careless spending."
    family_history = "Family history of mood disorders."
    diagnosis = "These symptoms, combined with the family history, strongly suggest a manic episode, characteristic of Bipolar Disorder."
    print(f"   - The patient presents with: {symptoms}")
    print(f"   - This clinical picture points towards a diagnosis of: {diagnosis}\n")

    # Step 2: Consider the treatment for the initial symptoms.
    print("Step 2: Identifying the likely treatment.")
    treatment_context = "The patient was prescribed a medication for these behavioral disturbances (mania)."
    likely_medication = "A first-line, common medication for Bipolar Disorder is Lithium, a mood stabilizer."
    print(f"   - {treatment_context}")
    print(f"   - Therefore, the prescribed medication is very likely to be Lithium.\n")

    # Step 3: Analyze the new symptom that developed after treatment.
    print("Step 3: Analyzing the subsequent symptom.")
    new_symptom = "Decreased interest in having sex (sexual dysfunction)."
    timing = "This symptom appeared *after* starting the new medication (likely Lithium)."
    print(f"   - The patient later developed: {new_symptom}")
    print(f"   - Crucially, the timing was: {timing}\n")

    # Step 4: Connect the treatment to the new symptom and evaluate the options.
    print("Step 4: Evaluating the causal link between the likely medication and the new symptom.")
    explanation = """
    We need to find a known side effect of Lithium that causes sexual dysfunction.
    - Lithium is well-known to cause Hypothyroidism (an underactive thyroid gland) in some patients.
    - Hypothyroidism is a well-known medical cause of sexual dysfunction, including decreased libido.
    - This creates a direct causal chain: Mania -> Lithium Treatment -> Hypothyroidism -> Sexual Dysfunction.
    - Other options involving heavy metals (Lead, Arsenic, etc.) do not explain this specific sequence of events, particularly the onset of sexual dysfunction *after* the new medication was started.
    """
    print(textwrap.indent(explanation, '   '))

    # Step 5: Final conclusion.
    print("Step 5: Conclusion.")
    final_answer = "A. Lithium induced hypothyroidism"
    print(f"   - The most comprehensive explanation that fits all the facts in the vignette is that the patient's sexual dysfunction was caused by hypothyroidism, which in turn was induced by the Lithium he was prescribed for his Bipolar Disorder.")
    print("\nFinal Answer Choice:")
    print(f"   {final_answer}")


if __name__ == "__main__":
    solve_clinical_case()
    print("<<<A>>>")