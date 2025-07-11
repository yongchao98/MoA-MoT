import textwrap

def explain_renal_decline_cause():
    """
    This script analyzes the provided clinical scenario to determine the lab parameter
    that would best indicate the cause of the patient's rapid renal decline.
    """
    
    print("Clinical Analysis Steps:")
    print("======================\n")
    
    # Step 1: Synthesize the clinical picture
    synthesis = (
        "1. Patient Presentation: The patient's history of facial rash, joint pain, fever, "
        "and pre-existing blood in the urine strongly suggests a systemic autoimmune disease, "
        "most likely Systemic Lupus Erythematosus (SLE)."
    )
    print(textwrap.fill(synthesis, 80))
    print("-" * 80)
    
    # Step 2: Identify the acute problem
    acute_event = (
        "2. Acute Deterioration: The rapid decline in kidney function—characterized by "
        "increased blood in the urine, decreased urine output, and headache—is a classic "
        "presentation of rapidly progressive glomerulonephritis (RPGN), a severe form of kidney inflammation."
    )
    print(textwrap.fill(acute_event, 80))
    print("-" * 80)

    # Step 3: Explain the underlying pathology
    pathology = (
        "3. Pathological Cause: In patients with SLE, this type of kidney failure is "
        "caused by the deposition of immune complexes in the kidney's glomeruli. These "
        "complexes are formed when autoantibodies, like anti-dsDNA, bind to antigens."
    )
    print(textwrap.fill(pathology, 80))
    print("-" * 80)

    # Step 4: Link pathology to a specific lab test
    lab_marker_explanation = (
        "4. Identifying the Best Lab Marker: The deposition of these immune complexes "
        "triggers the activation of the complement system. This process consumes complement "
        "proteins, primarily C3 and C4. Therefore, measuring their levels in the blood "
        "is a direct indicator of this pathological process. Low levels of C3 and C4 "
        "(hypocomplementemia) would be the most direct evidence of the immune-complex "
        "mediated damage that is causing the rapid decline in renal function."
    )
    print(textwrap.fill(lab_marker_explanation, 80))
    print("\n" + "=" * 80)

    # Conclusion
    conclusion = (
        "Conclusion: While markers like serum creatinine confirm kidney failure and anti-dsDNA "
        "antibodies point to high SLE activity, low C3 and C4 levels would most directly "
        "indicate the *cause* of the damage—the ongoing consumption of complement by immune "
        "complexes within the kidneys."
    )
    print(textwrap.fill(conclusion, 80))

explain_renal_decline_cause()