import textwrap

def explain_clinical_reasoning():
    """
    Analyzes the clinical scenario to identify the most indicative lab parameter.
    """

    print("Step 1: Analyzing the Patient's Chronic Symptoms")
    chronic_symptoms = "The patient's 7-year history of facial rash, joint pain, recurrent fever, and blood in the urine strongly suggests a systemic autoimmune disease, most likely Systemic Lupus Erythematosus (SLE)."
    print(textwrap.fill(chronic_symptoms, 80))
    print("-" * 80)

    print("Step 2: Understanding the Acute Deterioration")
    acute_event = "The rapid decline in kidney function after stopping corticosteroid treatment indicates a severe disease flare. This presentation is characteristic of rapidly progressive lupus nephritis, a life-threatening complication of SLE."
    print(textwrap.fill(acute_event, 80))
    print("-" * 80)

    print("Step 3: Identifying the Underlying Immunological Cause")
    mechanism = "The primary cause of kidney damage in lupus nephritis is the deposition of immune complexes (antibody-antigen molecules) in the kidney's filters (glomeruli). These complexes trigger a powerful inflammatory cascade by activating the complement system."
    print(textwrap.fill(mechanism, 80))
    print("-" * 80)

    print("Step 4: Connecting the Cause to a Specific Lab Test")
    conclusion = "When the complement system is activated so intensely, its protein components, particularly C3 and C4, are consumed faster than they can be produced. Therefore, a significant decrease in serum C3 and C4 levels is a direct and powerful indicator of this active, tissue-damaging immune process. This drop often precedes or occurs simultaneously with the most severe clinical signs of kidney failure."
    print(textwrap.fill(conclusion, 80))
    print("-" * 80)

    print("Final Conclusion:")
    final_answer = "While elevated anti-dsDNA antibodies confirm the lupus flare and rising creatinine confirms kidney failure, the parameter that best indicates the immunological *cause* of the rapid tissue destruction is decreased levels of complement proteins C3 and C4."
    print(textwrap.fill(final_answer, 80))

# Execute the function to provide the explanation.
explain_clinical_reasoning()