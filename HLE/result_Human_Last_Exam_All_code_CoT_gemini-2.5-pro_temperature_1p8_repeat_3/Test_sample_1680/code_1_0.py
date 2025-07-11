import operator

def analyze_patient_pathology():
    """
    Analyzes a clinical vignette to determine the best categorization for a patient's pathology.
    """

    case_facts = {
        "symptoms": ["memory loss", "forgets to feed himself", "weight loss", "disorientation"],
        "patient_statements": {
            "immediate_recall": "correct",
            "denial": True,
            "confabulation": "rare tapeworm for weight loss"
        },
        "history": {
            "positive": ["chronic venous insufficiency", "10 pack years smoking"],
            "negative": ["hypertension", "cirrhosis"]
        },
        "exam": "normal"
    }

    # Answer choices and their corresponding evidence mapping (1 for support, -1 for contradiction)
    analysis = {
        'A. Short-term memory': {
            'Memory loss is the primary complaint': 1,
            'Forgetting to eat and disorientation are signs of short-term memory deficits': 1,
            'These symptoms are the central feature of the clinical picture': 1
        },
        'B. Restrictive cardiomyopathy': {
            'Physical exam is normal, making a significant cardiac issue unlikely': -1,
            'No mention of typical symptoms like shortness of breath or edema (venous insufficiency is a distractor)': 0
        },
        'C. Hepatic encephalopathy': {
            'History explicitly states pertinent negative of cirrhosis, which is the primary cause': -1
        },
        'D. Parasitic infection': {
            'Claim made by a cognitively impaired patient, likely a confabulation to explain weight loss': -1,
            'No objective evidence presented': -1
        },
        'E. ATP depletion': {
            'This is a cellular-level mechanism, not a clinical category or syndrome': -1,
            'It is too general and non-specific to be a useful diagnosis': -1
        }
    }

    results = {}
    print("Analyzing each option based on the clinical vignette:\n")

    for option, evidence in analysis.items():
        score = sum(evidence.values())
        results[option] = score
        
        # Build the equation string
        equation_parts = [f"{val} ({reason})" for reason, val in evidence.items()]
        equation_str = " + ".join(equation_parts).replace("+ -", "- ")
        
        print(f"Analysis for {option}:")
        print(f"Final Equation: {equation_str} = {score}\n")

    # Determine the best option
    best_option = max(results.items(), key=operator.itemgetter(1))[0]
    
    print("--- Conclusion ---")
    print(f"The analysis shows that '{best_option}' has the highest relevance score.")
    print("The patient's core issue is a significant cognitive decline, with the most prominent feature being short-term memory impairment leading to functional decline. Other options are either contradicted by the case details or are not appropriate clinical categorizations.")

analyze_patient_pathology()