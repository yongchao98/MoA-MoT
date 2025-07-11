def get_lab_indicator():
    """
    Analyzes the clinical scenario to determine the most likely lab parameter
    indicating the cause of rapid renal decline.
    """
    # Patient history suggests Systemic Lupus Erythematosus (SLE).
    # The rapid deterioration into end-stage kidney disease points to a severe flare of lupus nephritis.
    # Lupus nephritis is caused by immune complexes depositing in the kidneys.
    # This process activates and consumes complement proteins.
    # Therefore, a decrease in complement levels is a direct marker of the disease process causing the kidney damage.
    
    best_indicator = "Decreased serum complement levels (C3 and C4)"
    explanation = (
        "Given the patient's history, which is highly suggestive of Systemic Lupus Erythematosus (SLE), "
        "the rapid decline in kidney function is characteristic of a severe lupus nephritis flare. "
        "The underlying cause of kidney damage in lupus nephritis is the deposition of immune complexes "
        "in the glomeruli. This process activates the complement cascade, leading to the consumption of "
        "complement proteins. Therefore, a significant decrease in serum complement levels (C3 and C4) "
        "would have been the best laboratory indicator reflecting the acute, immune-mediated "
        "pathology that caused the rapid renal failure."
    )
    
    print(explanation)

get_lab_indicator()