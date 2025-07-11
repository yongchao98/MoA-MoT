import textwrap

def find_lab_marker():
    """
    Analyzes the clinical scenario to determine the most indicative lab parameter
    for the cause of the patient's rapid renal decline.
    """

    # Clinical findings extracted from the scenario
    chronic_symptoms = ["facial rash", "joint pain", "recurrent fever", "hematuria"]
    acute_kidney_symptoms = ["significant increase in hematuria", "decreased urine output", "end-stage kidney disease"]
    trigger = "Inflammatory rebound after discontinuing corticosteroids"

    # Step-by-step reasoning
    print("Step 1: Assessing the Underlying Condition")
    print(textwrap.fill(
        f"The patient's long-standing history including {', '.join(chronic_symptoms)} is highly characteristic of the autoimmune disease Systemic Lupus Erythematosus (SLE).",
        80
    ))
    print("-" * 80)

    print("Step 2: Assessing the Acute Kidney Injury")
    print(textwrap.fill(
        f"The rapid deterioration of kidney function following an inflammatory rebound ({trigger}) points to a severe flare of SLE affecting the kidneys. This condition is known as Lupus Nephritis, which has progressed rapidly.",
        80
    ))
    print("-" * 80)

    print("Step 3: Identifying the Causative Lab Marker")
    print(textwrap.fill(
        "The question asks for the lab parameter indicating the CAUSE of the renal decline. In a lupus nephritis flare, the pathogenic mechanism involves the formation and deposition of immune complexes in the kidneys. These complexes are primarily formed by anti-double-stranded DNA (anti-dsDNA) antibodies.",
        80
    ))
    print("\n" + textwrap.fill(
        "Therefore, a high titer of anti-dsDNA antibodies is not just a marker of disease activity but is considered a key component of the immunological attack causing the kidney damage. While other markers like low complement (C3/C4) and high creatinine are expected, anti-dsDNA antibodies best indicate the direct cause of the immune-mediated assault on the kidneys.",
        80
    ))
    print("-" * 80)

    final_answer = "A high titer of anti-double-stranded DNA (anti-dsDNA) antibodies"
    print(f"Conclusion: The lab parameter that could have best indicated the cause of the rapid renal function decline is:\n{final_answer}")

if __name__ == "__main__":
    find_lab_marker()