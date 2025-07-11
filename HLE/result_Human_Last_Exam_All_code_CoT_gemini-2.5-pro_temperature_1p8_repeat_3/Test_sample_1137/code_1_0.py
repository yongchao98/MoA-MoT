import sys
# Redirect print to a string to control the final output format precisely.
# Note: The actual code execution in a user's environment will print to the console.
original_stdout = sys.stdout 
from io import StringIO
sys.stdout = captured_output = StringIO()

def solve_litigation_forum():
    """
    Analyzes and scores potential litigation forums based on a legal scenario.
    """
    # Criteria for evaluation derived from the problem description:
    # 1. Is it a court of first instance? (Originating Court)
    # 2. Does it have jurisdiction over the provincial commercial matter? (Jurisdiction)
    # 3. Is it suitable for a high-value, complex case? (Complexity)
    # 4. Does it meet the specific requirement for a speedy resolution? (Speed)

    forums = {
        "A. Ontario Court of Appeal":   {"Originating": 0, "Jurisdiction": 1, "Complexity": 1, "Speed": 0},
        "B. Commercial List":           {"Originating": 1, "Jurisdiction": 1, "Complexity": 2, "Speed": 2},
        "C. Superior Court of Justice": {"Originating": 1, "Jurisdiction": 1, "Complexity": 1, "Speed": 1},
        "D. Small Claims Court":        {"Originating": 1, "Jurisdiction": 1, "Complexity": 0, "Speed": 1},
        "E. Federal Court of Canada":   {"Originating": 1, "Jurisdiction": 0, "Complexity": 1, "Speed": 1}
    }

    # Assigning a higher score (2) for an ideal fit on Complexity and Speed
    
    best_forum = None
    max_score = -1
    results = {}

    print("--- Analysis of Litigation Forum Options ---")
    for forum, scores in forums.items():
        total_score = sum(scores.values())
        results[forum] = {"scores": scores, "total": total_score}
        if total_score > max_score:
            max_score = total_score
            best_forum = forum

    print("The dispute involves a complex, high-value commercial matter governed by provincial law, with a key need for a speedy resolution.")
    print("\nScoring based on suitability (0=unsuitable, 1=suitable, 2=ideal):")
    for forum, data in results.items():
        print(f"- {forum}: Score = {data['total']}")

    winning_data = results[best_forum]
    winning_scores = winning_data["scores"]

    print(f"\nThe best option is '{best_forum}' with a total score of {max_score}.")
    
    # Printing the final equation with each number as requested
    print("\nThe final equation for the winning choice is:")
    
    # We retrieve the individual score values for the printout.
    originating_score = winning_scores['Originating']
    jurisdiction_score = winning_scores['Jurisdiction']
    complexity_score = winning_scores['Complexity']
    speed_score = winning_scores['Speed']
    
    print(f"Total Score = Originating Court Score + Jurisdiction Score + Complexity Score + Speed Score")
    print(f"{max_score} = {originating_score} + {jurisdiction_score} + {complexity_score} + {speed_score}")
    
    # The final answer in the required format.
    print("\n<<<B>>>")


solve_litigation_forum()
# Reset stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())