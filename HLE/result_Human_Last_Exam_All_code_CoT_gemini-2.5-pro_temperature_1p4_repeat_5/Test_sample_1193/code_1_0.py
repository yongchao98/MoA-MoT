import operator

def diagnose_hypoxemia():
    """
    Analyzes a clinical scenario to determine the most likely cause of hypoxemia
    by scoring each possible diagnosis against the patient's data.
    """
    # 1. Clinical Data from the scenario
    time_post_op_days = 29
    oxygen_level_percent = 82
    procedure = "Whipple procedure"
    physical_findings = "Bilateral crackles and gasping for air"
    history = "Major surgery with blood transfusions"

    print(f"Analyzing patient case: {time_post_op_days} days post-{procedure}...")
    print(f"Key signs: Oxygen at {oxygen_level_percent}%, {physical_findings}.")
    print("-" * 30)

    # 2. Define potential diagnoses (Answer Choices)
    diagnoses = {
        'A': "Acute blood transfusion reaction",
        'B': "Iodine-related reaction",
        'C': "Sensitivity reaction",
        'D': "Sepsis",
        'E': "Myocyte necrosis",
        'F': "Respiratory deconditioning",
        'G': "Lung exhaustion",
        'H': "Air pollution sensitivity",
    }

    # 3. Score each diagnosis based on the clinical data
    scores = {key: 0 for key in diagnoses.keys()}
    reasoning = {key: [] for key in diagnoses.keys()}

    # Score based on timeframe (29 days)
    scores['A'] -= 2  # Acute reactions happen in hours, not 29 days.
    reasoning['A'].append(f"Timeframe mismatch (-2): Acute reaction is unlikely {time_post_op_days} days later.")
    scores['D'] += 1  # Post-op infection can develop and cause sepsis over weeks.
    reasoning['D'].append(f"Timeframe plausible (+1): Post-op infection leading to sepsis is a known complication within this {time_post_op_days}-day window.")

    # Score based on symptoms (Severe hypoxemia, crackles -> ARDS picture)
    scores['D'] += 2  # Sepsis is a classic cause of ARDS.
    reasoning['D'].append(f"Strongly explains symptoms (+2): Sepsis is a primary cause of Acute Respiratory Distress Syndrome (ARDS), matching the severe hypoxemia and crackles.")
    scores['E'] += 1  # Possible, but sepsis is a better fit post-op.
    reasoning['E'].append("Partially explains symptoms (+1): Can cause fluid in lungs, but sepsis/ARDS is a more direct fit for the overall context.")
    scores['F'] -= 2  # Does not explain acute severity or crackles.
    reasoning['F'].append("Does not explain symptoms (-2): Deconditioning causes weakness, not severe hypoxemia with crackles.")
    scores['G'] -= 2  # Not a recognized medical term for this condition.
    reasoning['G'].append("Invalid diagnosis (-2): 'Lung exhaustion' is not a standard medical term.")

    # Score based on context (Post-Whipple)
    scores['D'] += 2 # Whipple procedure has a high risk of infection.
    reasoning['D'].append(f"Matches context perfectly (+2): The {procedure} carries a high risk of abdominal infection, which is a common source of sepsis.")

    # 4. Determine the best diagnosis
    best_choice = max(scores.items(), key=operator.itemgetter(1))[0]

    print("Evaluation Results (Score):")
    for key, score in scores.items():
        print(f"  {key}: {diagnoses[key]:<35} | Score: {score}")
    print("-" * 30)

    # 5. Show the equation for the winning answer
    print(f"Final Analysis for the Highest Scoring Option: '{diagnoses[best_choice]}'")
    
    # Manually define the components of the score for the printout
    timeframe_score = 1
    symptoms_score = 2
    context_score = 2
    total_score = timeframe_score + symptoms_score + context_score

    print("\nReasoning:")
    for r in reasoning[best_choice]:
          print(f"- {r}")

    print("\nFinal Score Calculation:")
    print(f"Score from timeframe ({timeframe_score}) + Score from symptoms ({symptoms_score}) + Score from context ({context_score}) = {total_score}")


diagnose_hypoxemia()
<<<D>>>