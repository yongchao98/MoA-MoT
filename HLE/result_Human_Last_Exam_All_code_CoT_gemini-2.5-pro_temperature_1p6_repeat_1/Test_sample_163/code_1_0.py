import collections

def evaluate_surveillance_plans():
    """
    This function evaluates post-procedure surveillance plans for SFA stenting
    based on established medical guidelines.
    """

    # Patient Case: 65-year-old male, post-SFA stenting for 90% stenosis.
    # Goal: Select the most appropriate surveillance program.

    # Guidelines require two key components for optimal surveillance:
    # 1. Imaging Modality: Arterial Duplex is superior to ABI alone for detecting in-stent restenosis.
    # 2. Frequency: Surveillance should be frequent in the first year (e.g., 3, 6, 12 months).

    plans = {
        'A': {
            "description": "Regular visits with assessment for interval change in symptoms, vascular examination, and ABI measurement beginning in the immediate postprocedure period and at intervals for at least 2 years",
            "imaging": ["ABI"],
            "schedule": "Vague intervals for 2 years"
        },
        'B': {
            "description": "Regular visits with assessment for interval change in symptoms, vascular examination, and arterial duplex at 1 month, 3 months, and at month 12",
            "imaging": ["Arterial Duplex"],
            "schedule": "1, 3, 12 months"
        },
        'C': {
            "description": "Regular visits with assessment for interval change in symptoms, vascular examination, and ABI measurement at 3 months, 6 months, 9 months, and at month 12",
            "imaging": ["ABI"],
            "schedule": "3, 6, 9, 12 months"
        },
        'D': {
            "description": "Regular visits with assessment for interval change in symptoms, vascular examination, and arterial duplex at 3 months, 6 months, 12 months, and 2 years",
            "imaging": ["Arterial Duplex"],
            "schedule": "3, 6, 12 months, and 2 years"
        },
        'E': {
            "description": "Annual visits with assessment for interval change in symptoms, vascular examination, ABI measurement, and arterial duplex",
            "imaging": ["ABI", "Arterial Duplex"],
            "schedule": "Annual"
        }
    }

    best_option = None
    max_score = -1

    print("Evaluating surveillance options based on medical best practices:\n")

    for option, details in plans.items():
        score = 0
        reasons = []

        # Criterion 1: Check for Duplex Ultrasound
        if "Arterial Duplex" in details["imaging"]:
            score += 2
            reasons.append("Includes Duplex Ultrasound (Correct Modality).")
        else:
            reasons.append("Fails to include Duplex Ultrasound, relying only on ABI.")

        # Criterion 2: Check for appropriate frequency
        if "6 months" in details["schedule"] and "12 months" in details["schedule"]:
            score += 1
            reasons.append("Includes standard 6 and 12-month follow-ups (Correct Frequency).")
        elif "Annual" in details["schedule"]:
             reasons.append("Frequency is too low for the first year.")
        else:
            reasons.append("Frequency is not ideal (e.g., missing a 6-month check or too vague).")
            
        print(f"Option {option}:")
        print(f"  - Plan: {details['description']}")
        print(f"  - Analysis: {' '.join(reasons)}")
        print(f"  - Score: {score}\n")

        if score > max_score:
            max_score = score
            best_option = option

    print("--- Conclusion ---")
    print(f"The best option is '{best_option}' because it combines the most appropriate imaging modality (Arterial Duplex) with the standard surveillance frequency (3, 6, and 12 months) to effectively monitor for in-stent restenosis during the critical first year.")
    
evaluate_surveillance_plans()
<<<D>>>