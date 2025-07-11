import sys

# Define a function to evaluate the diagnoses based on clinical findings
def evaluate_diagnosis():
    """
    This function uses a scoring system to determine the most likely diagnosis
    based on the key findings from the case study.
    """

    # Initialize scores for each diagnosis
    # A. Colonic perforation, B. Lower GI bleeding, C. Splenic laceration, D. Postpolypectomy syndrome
    scores = {
        'A': 0,
        'B': 0,
        'C': 0,
        'D': 0
    }

    diagnoses = {
        'A': "Colonic perforation",
        'B': "Lower GI bleeding",
        'C': "Splenic laceration",
        'D': "Postpolypectomy syndrome"
    }

    # Store the components of the final winning equation
    equation_components = []

    # 1. Pain Location: LUQ pain and left-shoulder pain (Kehr's sign)
    # This strongly suggests splenic injury irritating the diaphragm.
    scores['C'] += 3  # Classic presentation
    scores['A'] += 1  # Possible but less specific
    scores['B'] += 0  # Wrong location
    scores['D'] += 0  # Unlikely
    equation_components.append("3 (for LUQ and shoulder pain)")

    # 2. Hemodynamic Instability: Rapid, massive drop in hemoglobin.
    # This points to acute, major hemorrhage, characteristic of a solid organ injury.
    scores['C'] += 3  # Classic presentation
    scores['B'] += 1  # Possible but less common to be this massive without polypectomy
    scores['A'] += 1  # Possible but often presents more with infection signs
    scores['D'] += 0  # Not a bleeding syndrome
    equation_components.append("3 (for rapid hemoglobin drop)")

    # 3. Procedure Detail: No polypectomy was performed.
    # This fact rules out postpolypectomy syndrome entirely.
    scores['D'] -= 10 # This is a rule-out criterion.
    # It makes simple lower GI bleeding less likely as well.
    scores['B'] -= 1
    scores['A'] += 1  # Not dependent on polypectomy
    scores['C'] += 1  # Not dependent on polypectomy
    equation_components.append("1 (for occurring without polypectomy)")
    
    # 4. Procedure Detail: "Difficult" colonoscopy
    # This is a known risk factor for splenic injury and perforation due to torque.
    scores['C'] += 2  # Strong risk factor
    scores['A'] += 2  # Strong risk factor
    scores['B'] += 0  # Not a direct risk factor
    scores['D'] += 0  # Irrelevant
    equation_components.append("2 (for occurring after a difficult procedure)")


    # Find the diagnosis with the highest score
    most_likely_diagnosis_code = max(scores, key=scores.get)
    highest_score = scores[most_likely_diagnosis_code]
    most_likely_diagnosis_name = diagnoses[most_likely_diagnosis_code]

    # Print the "equation" for the winning diagnosis as requested
    equation_string = " + ".join(equation_components)
    
    print(f"Analysis complete. The most likely diagnosis is {most_likely_diagnosis_name}.")
    print("\nScoring Equation:")
    print(f"{most_likely_diagnosis_name} (C) Score = {equation_string} = {highest_score}")
    print("\nExplanation:")
    print("The patient's presentation of left upper quadrant pain, referred left shoulder pain (Kehr's sign), and hemodynamic shock from a massive drop in hemoglobin after a difficult colonoscopy is the classic triad for splenic laceration.")
    print("\nFinal Answer Code:")
    sys.stdout.write("<<<C>>>")

# Run the evaluation
evaluate_diagnosis()