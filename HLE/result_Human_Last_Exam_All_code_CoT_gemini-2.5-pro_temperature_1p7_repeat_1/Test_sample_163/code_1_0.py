import sys

def solve_medical_surveillance():
    """
    This script evaluates the best post-procedure surveillance program for a patient
    after superficial femoral artery (SFA) stenting by scoring each option based on
    current clinical best practices.
    """

    options = {
        'A': {'uses_duplex': False, 'timing': 'non-standard', 'text': 'Regular visits... ABI measurement beginning in the immediate postprocedure period and at intervals for at least 2 years'},
        'B': {'uses_duplex': True, 'timing': 'suboptimal', 'text': 'Regular visits... arterial duplex at 1 month, 3 months, and at month 12'},
        'C': {'uses_duplex': False, 'timing': 'suboptimal', 'text': 'Regular visits... ABI measurement at 3 months, 6 months, 9 months, and at month 12'},
        'D': {'uses_duplex': True, 'timing': 'optimal', 'text': 'Regular visits... arterial duplex at 3 months, 6 months, 12 months, and 2 years'},
        'E': {'uses_duplex': True, 'timing': 'infrequent', 'text': 'Annual visits... ABI measurement, and arterial duplex'}
    }

    scores = {}
    best_option = ''
    max_score = -sys.maxsize

    print("Evaluating surveillance options based on two key criteria:")
    print("1. Use of Arterial Duplex Ultrasound (the gold standard for detecting in-stent restenosis).")
    print("2. Appropriate Surveillance Timing (frequent checks in the first year).\n")

    # Define score components
    duplex_score_component = 0
    timing_score_component = 0

    for key, details in options.items():
        # Criterion 1: Use of Duplex Ultrasound
        duplex_score = 10 if details['uses_duplex'] else -20

        # Criterion 2: Timing of Surveillance
        if details['timing'] == 'optimal':
            timing_score = 10  # Standard of care: 3, 6, 12 months is key
        elif details['timing'] == 'infrequent':
            timing_score = -5 # Annual is too infrequent for the first year
        elif details['timing'] == 'suboptimal':
            timing_score = 2   # Better than nothing but misses key timepoints
        else:
            timing_score = 0

        total_score = duplex_score + timing_score
        scores[key] = total_score

        if total_score > max_score:
            max_score = total_score
            best_option = key
            duplex_score_component = duplex_score
            timing_score_component = timing_score

    print("--- Scoring Results ---")
    for key in options.keys():
        print(f"Option {key}: Score = {scores[key]}")

    print("\n--- Conclusion ---")
    print(f"The optimal choice is Option '{best_option}' with the highest score of {max_score}.")
    print("This choice correctly includes arterial duplex at standard, frequent intervals during the first year (3, 6, 12 months), which is critical for detecting early stent failure.\n")
    print("The final score calculation for the best option is:")
    print(f"Total Score = (Duplex Use Score) + (Timing Score)")
    print(f"{max_score} = ({duplex_score_component}) + ({timing_score_component})\n")
    
    print("Therefore, the most appropriate surveillance program is:")
    print(f"D. {options['D']['text']}")

solve_medical_surveillance()