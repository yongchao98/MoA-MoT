import pandas as pd

def find_best_litigation_forum():
    """
    Analyzes a legal scenario to determine the best litigation forum.

    This function models the decision-making process by scoring potential courts
    based on criteria derived from the case facts.
    """

    # 1. Define the requirements of the legal case
    case_requirements = {
        'needs_trial_court': True,
        'is_high_value': True,
        'is_complex_commercial': True,
        'needs_speed': True,
        'is_provincial_jurisdiction': True
    }

    # 2. Define the attributes of each potential forum
    forums = {
        'A. Ontario Court of Appeal': {
            'is_trial_court': False, 'has_jurisdiction': True, 'is_expedited': False, 'has_commercial_expertise': True
        },
        'B. Commercial List': {
            'is_trial_court': True, 'has_jurisdiction': True, 'is_expedited': True, 'has_commercial_expertise': True
        },
        'C. Superior Court of Justice': {
            'is_trial_court': True, 'has_jurisdiction': True, 'is_expedited': False, 'has_commercial_expertise': False
        },
        'D. Small Claims Court': {
            'is_trial_court': True, 'has_jurisdiction': False, 'is_expedited': False, 'has_commercial_expertise': False
        },
        'E. Federal Court of Canada': {
            'is_trial_court': True, 'has_jurisdiction': False, 'is_expedited': False, 'has_commercial_expertise': False
        }
    }
    
    # 3. Score each forum based on how well it matches the case requirements
    # We assign different weights to each criterion based on importance.
    # Speed and expertise are the most critical differentiators in this scenario.
    scores = {}
    score_equations = {}

    weights = {
        'trial_court_fit': 1,
        'jurisdiction_fit': 10, # A mismatch is a fatal flaw
        'expertise_fit': 3,     # Highly important for this complex case
        'speed_fit': 3          # A key stated goal for the client
    }
    
    for name, attrs in forums.items():
        score = 0
        equation_parts = []

        # A. Must be a trial court
        if attrs['is_trial_court'] == case_requirements['needs_trial_court']:
            score += weights['trial_court_fit']
            equation_parts.append(str(weights['trial_court_fit']))
        else:
            equation_parts.append("0")


        # B. Must have jurisdiction (monetary and subject matter)
        if attrs['has_jurisdiction'] == case_requirements['is_provincial_jurisdiction']:
            score += weights['jurisdiction_fit']
            equation_parts.append(str(weights['jurisdiction_fit']))
        else:
            # If jurisdiction is wrong, it's not a viable option.
            score = -1 
            equation_parts = ["-1 (Fatal Flaw: No Jurisdiction)"]

        # C. Specialized commercial expertise is best
        if attrs['has_commercial_expertise'] and case_requirements['is_complex_commercial'] and score > 0:
            score += weights['expertise_fit']
            equation_parts.append(str(weights['expertise_fit']))
        else:
            equation_parts.append("0")
            
        # D. Expedited process is required
        if attrs['is_expedited'] and case_requirements['needs_speed'] and score > 0:
            score += weights['speed_fit']
            equation_parts.append(str(weights['speed_fit']))
        else:
            equation_parts.append("0")
        
        scores[name] = score
        score_equations[name] = " + ".join(equation_parts)

    # 4. Determine the best option and print the results
    best_forum = max(scores, key=scores.get)
    
    print("Analysis of Litigation Forum Options:\n")
    print("The case involves a high-value, complex commercial dispute under Ontario law, with a key requirement for a speedy resolution.\n")
    
    df = pd.DataFrame.from_dict(forums, orient='index')
    df['Score'] = pd.Series(scores)
    print(df)
    
    print("\n-------------------------------------------------\n")
    print(f"Conclusion: The best available choice is '{best_forum}'.\n")
    print("Reasoning: This forum is a trial court with the correct jurisdiction. Crucially, it is specifically designed for complex commercial cases and offers an expedited process, directly addressing all of RE1's needs.\n")
    
    # Final equation as requested
    final_equation = score_equations[best_forum]
    final_score = scores[best_forum]
    print("Final Calculation for the Best Option:")
    print(f"Score for '{best_forum}' is derived from fitting the criteria (Trial Court + Jurisdiction + Expertise + Speed).")
    print(f"Equation: {final_equation} = {final_score}")

find_best_litigation_forum()