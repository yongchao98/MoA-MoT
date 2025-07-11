import textwrap

def solve_diagnostic_case():
    """
    Analyzes a clinical case to determine the best next diagnostic step.
    """

    # Step 1: Analyze the patient's presentation from the provided text.
    case_summary = {
        "Complaint": "Itchy rash in both axillae",
        "Key Finding": "Rash on the posterior border of axillary folds, SPARING the axillary vaults.",
        "History Clue": "Started a new workout program 3 weeks ago (implying new workout clothes).",
        "Suspected Diagnosis": "Allergic contact dermatitis, specifically textile dermatitis from clothing dyes or resins."
    }

    # Step 2: Define and evaluate the diagnostic options.
    options = {
        "A": "Skin biopsy",
        "B": "KOH preparation",
        "C": "Topical steroid",
        "D": "Patch test",
        "E": "None of the above"
    }

    reasoning = {
        "A": "A skin biopsy is invasive and typically reserved for when the diagnosis is uncertain or to rule out more serious conditions. The clinical picture here strongly points to contact dermatitis, making a biopsy unnecessary at this stage.",
        "B": "A KOH prep tests for fungal infections. While fungus can cause rashes in the armpits, the specific distribution (sparing the vault) is classic for textile dermatitis, not a fungal infection.",
        "C": "A topical steroid is a TREATMENT for inflammation and itching, not a DIAGNOSTIC step. The question asks for the next step in diagnosis.",
        "D": "A patch test is the gold standard for identifying the specific causative allergen in allergic contact dermatitis. Since the suspected cause is a chemical (dye or resin) in the patient's new workout clothes, this test would confirm the diagnosis and identify the trigger.",
        "E": "Option D is a valid and appropriate next step."
    }

    # Step 3: Print the analysis and conclusion.
    print("Analyzing the Clinical Case...")
    print("-" * 30)
    for key, value in case_summary.items():
        print(f"{key}: {value}")

    print("\nEvaluating the Answer Choices...")
    print("-" * 30)
    for option, description in options.items():
        if option != "E":
            print(f"Choice {option} ({description}):")
            # Using textwrap for better formatting of long reasoning strings
            wrapped_text = textwrap.fill(reasoning[option], width=70, initial_indent="  ", subsequent_indent="  ")
            print(wrapped_text)
            print("-" * 10)

    print("\nConclusion:")
    print("The patient's history and the physical exam strongly suggest allergic contact dermatitis from clothing (textile dermatitis). The most definitive step to confirm this diagnosis and identify the specific allergen is a patch test.")

    # Final Answer
    final_answer = "D"
    print(f"\nTherefore, the best next step in diagnosis is the Patch Test.")
    print(f"\n<<<D>>>")

solve_diagnostic_case()