import sys

def analyze_litigation_forum():
    """
    Analyzes the best litigation forum for a given legal scenario by scoring each option.
    """
    # 1. Define the key characteristics of the case from the problem description.
    case_factors = {
        'is_initial_claim': 1,      # The case is starting, not an appeal.
        'is_provincial_jurisdiction': 1, # The dispute is between Ontario companies over Ontario property.
        'is_complex_commercial': 1,  # Involves joint ventures, undertakings, and financial disputes.
        'is_high_value': 1,         # Six large commercial properties are at stake.
        'desires_speed': 1          # RE1 wants the claim resolved in the shortest time.
    }

    # 2. Define the attributes of each potential forum.
    # Scores are binary (1=Yes, 0=No) for matching attributes.
    forums = {
        "A. Ontario Court of Appeal": {
            'handles_initial_claim': 0, 'is_provincial': 1, 'handles_complex': 1,
            'handles_high_value': 1, 'is_expedited': 0, 'is_specialized': 0
        },
        "B. Commercial List": {
            'handles_initial_claim': 1, 'is_provincial': 1, 'handles_complex': 1,
            'handles_high_value': 1, 'is_expedited': 1, 'is_specialized': 1 # Specialized in commercial matters
        },
        "C. Superior Court of Justice": {
            'handles_initial_claim': 1, 'is_provincial': 1, 'handles_complex': 1,
            'handles_high_value': 1, 'is_expedited': 0, 'is_specialized': 0
        },
        "D. Small Claims Court": {
            'handles_initial_claim': 1, 'is_provincial': 1, 'handles_complex': 0,
            'handles_high_value': 0, 'is_expedited': 1, 'is_specialized': 0
        },
        "E. Federal Court of Canada": {
            'handles_initial_claim': 1, 'is_provincial': 0, 'handles_complex': 1,
            'handles_high_value': 1, 'is_expedited': 0, 'is_specialized': 0
        }
    }

    best_forum = ""
    max_score = -1
    best_forum_equation = ""

    # 3. Calculate a score for each forum.
    print("Analysis of Litigation Forums:")
    for name, attributes in forums.items():
        score = 0
        equation_parts = []

        # Calculate score based on matching factors
        score_initial = attributes['handles_initial_claim'] * case_factors['is_initial_claim']
        score_jurisdiction = attributes['is_provincial'] * case_factors['is_provincial_jurisdiction']
        score_complexity = attributes['handles_complex'] * case_factors['is_complex_commercial']
        score_value = attributes['handles_high_value'] * case_factors['is_high_value']
        score_speed = attributes['is_expedited'] * case_factors['desires_speed']
        score_specialized = attributes['is_specialized'] # Bonus for specialization

        score = score_initial + score_jurisdiction + score_complexity + score_value + score_speed + score_specialized
        
        current_equation = (f"Score = {score_initial} (initial) + {score_jurisdiction} (jurisdiction) + "
                            f"{score_complexity} (complex) + {score_value} (value) + "
                            f"{score_speed} (speed) + {score_specialized} (specialized) = {score}")

        print(f"- {name}: {current_equation}")

        if score > max_score:
            max_score = score
            best_forum = name
            best_forum_equation_values = [
                score_initial, score_jurisdiction, score_complexity,
                score_value, score_speed, score_specialized, score
            ]

    # 4. Print the final recommendation and its scoring equation.
    print("\n--- Conclusion ---")
    print(f"The best choice is '{best_forum}' with a score of {max_score}.")
    print("\nReasoning: The case is a complex, high-value commercial dispute within Ontario provincial jurisdiction, and the client desires a speedy resolution. The Commercial List is a specialized section of the Superior Court designed specifically for such cases, offering expert judges and active case management to ensure efficiency.")
    
    # Print the final equation with each number as requested
    final_values = best_forum_equation_values
    print("\nFinal Recommended Forum Equation:")
    print(f"'{best_forum}' Score = {final_values[0]} + {final_values[1]} + {final_values[2]} + {final_values[3]} + {final_values[4]} + {final_values[5]} = {final_values[6]}")

if __name__ == '__main__':
    analyze_litigation_forum()