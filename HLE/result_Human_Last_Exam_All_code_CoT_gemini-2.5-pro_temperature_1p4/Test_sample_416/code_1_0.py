def solve_diagnosis():
    """
    This function models the differential diagnosis process for the given clinical case
    by scoring potential diagnoses based on key findings.
    """
    # Potential diagnoses with initial scores
    scores = {
        "Osteoarthritis": 0,
        "Charcot Arthropathy": 0,
        "Septic Arthritis": 0,
        "Chronic osteomyelitis": 0,
        "Pseudogout": 0,
    }

    # Points assigned for each finding. Positive points support the diagnosis, negative points contradict it.
    # The structure is: finding_name: {diagnosis: points, ...}
    findings_points = {
        "Acute red, swollen joint": {
            "Osteoarthritis": -2, "Charcot Arthropathy": 5, "Septic Arthritis": 5, "Chronic osteomyelitis": 1, "Pseudogout": 5
        },
        "Triggered by minor trauma (long walk)": {
            "Osteoarthritis": 2, "Charcot Arthropathy": 5, "Septic Arthritis": 0, "Chronic osteomyelitis": 0, "Pseudogout": 1
        },
        "Initial X-rays negative": {
            "Osteoarthritis": -1, "Charcot Arthropathy": 5, "Septic Arthritis": 2, "Chronic osteomyelitis": -5, "Pseudogout": 2
        },
        "Failed treatment with NSAIDs/steroids": {
            "Osteoarthritis": -3, "Charcot Arthropathy": 5, "Septic Arthritis": -10, "Chronic osteomyelitis": -2, "Pseudogout": -10
        },
        "Synovial fluid: NO organisms or WBCs": {
            "Osteoarthritis": 2, "Charcot Arthropathy": 15, "Septic Arthritis": -20, "Chronic osteomyelitis": 2, "Pseudogout": 2
        },
        "Synovial fluid: NO crystals": {
            "Osteoarthritis": 1, "Charcot Arthropathy": 5, "Septic Arthritis": 1, "Chronic osteomyelitis": 1, "Pseudogout": -20
        },
    }

    # Keep track of the points added for the final "equation" output
    calculation_steps = {diag: [] for diag in scores}

    # Calculate scores by iterating through findings
    for finding, points in findings_points.items():
        for diagnosis, point_value in points.items():
            if diagnosis in scores:
                scores[diagnosis] += point_value
                calculation_steps[diagnosis].append(str(point_value))

    # Determine the most likely diagnosis
    most_likely_diagnosis = max(scores, key=scores.get)

    print("Differential Diagnosis Score Calculation:\n")
    for diagnosis, steps in calculation_steps.items():
        # The equation shows the points added or subtracted for each key clinical finding
        equation_str = " + ".join(steps).replace("+ -", "- ")
        final_score = scores[diagnosis]
        print(f"{diagnosis}:")
        print(f"Final Score = {final_score}")
        print(f"Calculation: Base(0) + {equation_str} = {final_score}\n")

    print(f"Based on the scoring model, the most likely diagnosis is: {most_likely_diagnosis}")

solve_diagnosis()