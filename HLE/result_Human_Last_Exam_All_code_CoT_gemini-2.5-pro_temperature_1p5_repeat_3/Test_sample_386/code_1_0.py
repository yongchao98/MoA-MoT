import pandas as pd

def evaluate_culture_protocol(options):
    """
    Evaluates cell culture protocols based on standard biological laboratory practices.

    Args:
        options (dict): A dictionary where keys are option letters and values are
                        dictionaries containing protocol parameters.

    Returns:
        str: The letter of the most plausible option.
    """
    scores = {}
    explanations = {}

    for option, params in options.items():
        score = 0
        explanation = []

        # 1. Evaluate Serum Concentration
        serum_concentration = params.get("serum_concentration", 0)
        is_serum_free = params.get("is_serum_free", False)

        if is_serum_free and serum_concentration > 0:
            score -= 2
            explanation.append(f"FAIL: '{serum_concentration}% serum-free medium' is a contradictory and invalid term.")
        elif not is_serum_free and 8 <= serum_concentration <= 15:
            score += 1
            explanation.append(f"PASS: {serum_concentration}% FBS is a standard concentration for fibroblast culture.")
        else:
            score -= 1
            explanation.append(f"FAIL: {serum_concentration}% FBS is outside the typical range (8-15%) or improperly defined.")

        # 2. Evaluate Antibiotic Concentration
        antibiotic_concentration = params.get("antibiotic_concentration", 0)
        if antibiotic_concentration == 1:
            score += 1
            explanation.append(f"PASS: {antibiotic_concentration}% antibiotics is the standard concentration.")
        else:
            score -= 1
            explanation.append(f"FAIL: {antibiotic_concentration}% antibiotics is not standard; 5% is likely toxic.")

        # 3. Evaluate Cell Adherence
        adherence = params.get("adherence")
        if adherence is True:
            score += 1
            explanation.append("PASS: Adherence to the flask is expected for a successful fibroblast culture.")
        elif adherence is False:
            score -= 2
            explanation.append("FAIL: Non-adherence indicates cell death or culture failure.")

        # 4. Evaluate Biological Premise
        premise = params.get("premise_ok", True)
        if not premise:
            score -= 2
            explanation.append("FAIL: The biological premise (e.g., cell source, differentiation pathway) is incorrect.")

        scores[option] = score
        explanations[option] = explanation

    # --- Analysis Output ---
    print("--- Analysis of Answer Choices ---")
    for option in options.keys():
        print(f"\nEvaluating Choice {option}:")
        for line in explanations[option]:
            print(f"- {line}")
        print(f"Total Plausibility Score: {scores[option]}")

    best_option = max(scores, key=scores.get)
    print("\n--- Conclusion ---")
    print(f"Choice {best_option} is the most biologically and procedurally correct statement.")

    return best_option

# Define the parameters from each answer choice
choices = {
    "A": {
        "description": "Stromal cells... prevented themselves from adhering... with 12% FBS and 1% antibiotics.",
        "serum_concentration": 12,
        "is_serum_free": False,
        "antibiotic_concentration": 1,
        "adherence": False
    },
    "B": {
        "description": "Limbal cell explants in... 10% serum-free medium and 5% antibiotic...",
        "serum_concentration": 10,
        "is_serum_free": True,
        "antibiotic_concentration": 5,
        "adherence": None,
        "premise_ok": False # Limbal cells are incorrect source for fibroblasts
    },
    "C": {
        "description": "Stromal cells to myofibroblasts in... 10% of the FBS and 1% antibiotic, and they adhered...",
        "serum_concentration": 10,
        "is_serum_free": False,
        "antibiotic_concentration": 1,
        "adherence": True
    },
    "D": {
        "description": "Limbal cells to proliferate into myofibroblasts...",
        "serum_concentration": 0, # Not specified, assuming none for scoring
        "is_serum_free": False,
        "antibiotic_concentration": 1,
        "adherence": None,
        "premise_ok": False # Incorrect differentiation pathway
    },
    "E": {
        "description": "Stromal cells... using 11% serum-free medium and 1% antibiotics...",
        "serum_concentration": 11,
        "is_serum_free": True,
        "antibiotic_concentration": 1,
        "adherence": None
    }
}

# Run the evaluation and find the best choice
correct_answer = evaluate_culture_protocol(choices)

# As requested, output the numbers from the final correct equation
final_serum = choices[correct_answer]['serum_concentration']
final_antibiotic = choices[correct_answer]['antibiotic_concentration']

print("\nThe final equation for the culture medium components is:")
print(f"Final Medium = Base Medium + {final_serum}% FBS + {final_antibiotic}% antibiotic")

# Final answer in the required format
print(f"\n\n<<<C>>>")