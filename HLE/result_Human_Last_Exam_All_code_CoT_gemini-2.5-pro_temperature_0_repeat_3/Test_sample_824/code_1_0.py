def find_best_lab_indicator():
    """
    Analyzes the clinical scenario to determine the most indicative lab parameter.

    The patient's history and symptoms strongly suggest Systemic Lupus Erythematosus (SLE)
    with lupus nephritis. The rapid renal decline after stopping steroids points to a
    severe disease flare.

    The cause of kidney damage in lupus nephritis is the deposition of immune complexes,
    which activates and consumes complement proteins. Therefore, measuring complement
    levels is key.
    """

    # The most direct indicator of the immune-mediated cause of the renal decline.
    best_indicator = "Complement levels (C3 and C4)"

    # Explanation
    explanation = (
        "The patient's rapid renal decline is characteristic of a severe flare of lupus nephritis. "
        "The underlying cause of the damage is an intense immunological attack on the kidneys. "
        "This process involves the deposition of immune complexes which activate and consume proteins of the complement system. "
        "Therefore, a significant decrease in complement levels, specifically C3 and C4, would be the best laboratory indicator "
        "of the active disease process causing the rapid decline in kidney function."
    )

    print(explanation)
    print("\n---")
    print(f"The lab parameter that could have best indicated the cause of the rapid renal function decline is: {best_indicator}")

find_best_lab_indicator()