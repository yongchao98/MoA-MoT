import pandas as pd

def analyze_clinical_case():
    """
    Analyzes a clinical vignette to find the most likely diagnosis by scoring
    each possibility against the provided evidence.
    """
    # Key findings from the case vignette
    patient_profile = {
        'initial_symptoms': "Agitation, difficulty sleeping, hypersexuality, careless spending (classic manic symptoms).",
        'family_history': "Mood disorders, predisposing to Bipolar Disorder.",
        'likely_medication': "Symptoms warrant a mood stabilizer, with Lithium being a common choice.",
        'subsequent_symptom': "Decreased sex drive (sexual dysfunction) developing *after* starting the new medication.",
        'occupational_history': "Metal smelting (raises possibility of heavy metal toxicity, but may be a distractor)."
    }

    # Define the answer choices and a scoring rubric
    options = {
        'A': 'Lithium induced hypothyroidism',
        'B': 'Arsenic induced Renal Dysfunction',
        'C': 'Mercury induced Renal Dysfunction',
        'D': 'Lead induced Sexual dysfunction',
        'E': 'Manganese induced Renal Dysfunction'
    }

    # Scoring criteria based on explanatory power for the sequence of events
    # Each criterion explains a part of the causal chain.
    scoring_points = {
        'explains_mania': 1,
        'explains_prescription': 1,
        'explains_sexual_dysfunction': 1,
        'explains_timing_of_dysfunction': 2, # Critical for connecting medication to the side effect
        'fits_history': 1,
    }

    # Evaluate each option
    scores = {}
    equations = {}

    # --- Scoring Logic ---

    # Option A: Lithium -> Hypothyroidism -> Dysfunction
    s_a1 = scoring_points['explains_mania'] # Assumes Bipolar, which Lithium treats
    s_a2 = scoring_points['explains_prescription'] # It IS the implied prescription
    s_a3 = scoring_points['explains_sexual_dysfunction'] # Hypothyroidism causes low libido
    s_a4 = scoring_points['explains_timing_of_dysfunction'] # Side effect develops over time
    s_a5 = scoring_points['fits_history'] # Fits family hx of mood disorders
    scores['A'] = s_a1 + s_a2 + s_a3 + s_a4 + s_a5
    equations['A'] = f"{s_a1} (explains mania) + {s_a2} (explains prescription) + {s_a3} (explains dysfunction) + {s_a4} (explains timing) + {s_a5} (fits family hx)"

    # Option B: Arsenic
    s_b1, s_b2, s_b4 = 0, 0, 0 # Doesn't explain mania, prescription, or timing
    s_b3 = scoring_points['explains_sexual_dysfunction'] # Renal issues can cause dysfunction
    s_b5 = scoring_points['fits_history'] # Fits occupational history
    scores['B'] = s_b1 + s_b2 + s_b3 + s_b4 + s_b5
    equations['B'] = f"{s_b1} (explains mania) + {s_b2} (explains prescription) + {s_b3} (explains dysfunction) + {s_b4} (explains timing) + {s_b5} (fits job hx)"

    # Option C: Mercury
    s_c1, s_c2, s_c4 = 0, 0, 0 # Doesn't explain mania, prescription, or timing
    s_c3 = scoring_points['explains_sexual_dysfunction'] # Renal issues can cause dysfunction
    s_c5 = scoring_points['fits_history'] # Fits occupational history
    scores['C'] = s_c1 + s_c2 + s_c3 + s_c4 + s_c5
    equations['C'] = f"{s_c1} (explains mania) + {s_c2} (explains prescription) + {s_c3} (explains dysfunction) + {s_c4} (explains timing) + {s_c5} (fits job hx)"

    # Option D: Lead
    s_d1, s_d2, s_d4 = 0, 0, 0 # Contradicts initial hypersexuality and doesn't explain timing
    s_d3 = scoring_points['explains_sexual_dysfunction'] # Lead can cause dysfunction
    s_d5 = scoring_points['fits_history'] # Fits occupational history
    scores['D'] = s_d1 + s_d2 + s_d3 + s_d4 + s_d5
    equations['D'] = f"{s_d1} (explains mania) + {s_d2} (explains prescription) + {s_d3} (explains dysfunction) + {s_d4} (explains timing) + {s_d5} (fits job hx)"
    
    # Option E: Manganese
    s_e1, s_e2, s_e4 = 0, 0, 0 # "Manganese madness" is a poor fit; doesn't explain prescription or timing
    s_e3 = scoring_points['explains_sexual_dysfunction'] # Renal issues can cause dysfunction
    s_e5 = scoring_points['fits_history'] # Fits occupational history
    scores['E'] = s_e1 + s_e2 + s_e3 + s_e4 + s_e5
    equations['E'] = f"{s_e1} (explains mania) + {s_e2} (explains prescription) + {s_e3} (explains dysfunction) + {s_e4} (explains timing) + {s_e5} (fits job hx)"

    # --- Print the analysis and result ---
    print("--- Clinical Case Analysis ---")
    for key, value in patient_profile.items():
        print(f"{key.replace('_', ' ').title()}: {value}")

    print("\n--- Scoring Each Potential Cause ---")
    for option_key in options:
        full_option_name = options[option_key]
        equation_str = equations[option_key]
        total_score = scores[option_key]
        print(f"Score for '{option_key}. {full_option_name}':")
        # Final equation output for each option
        print(f"    Equation: {equation_str} = {total_score}")

    # Determine the best option
    best_option_key = max(scores, key=scores.get)
    best_option_name = options[best_option_key]

    print("\n--- Conclusion ---")
    print(f"The root cause that best explains the entire series of events is:")
    print(f"'{best_option_name}', with the highest score of {scores[best_option_key]}.")
    print("\nRationale: The patient's initial symptoms are classic for a manic episode, justifying treatment with a mood stabilizer like Lithium. A well-known side effect of Lithium is hypothyroidism, which in turn is a common cause of decreased libido. This explanation cohesively links the patient's history, initial symptoms, treatment, and subsequent development of sexual dysfunction.")

if __name__ == "__main__":
    analyze_clinical_case()
    print("\n<<<A>>>")
