import sys

def solve_surgical_dilemma():
    """
    Analyzes surgical options for a stuck laparoscopic stapler
    and determines the safest course of action.
    """
    options = {
        'A': {
            "description": "Insert a laparoscopic stapler and staple to the right of the appendix resecting the appendix along with some of the cecum.",
            "risk": 7,  # High risk: Unnecessarily resects healthy bowel (cecum).
            "effectiveness": 6 # Effective in removal, but overly aggressive.
        },
        'B': {
            "description": "Insert a laparoscopic instrument like a flat grasper and get it between the jaws of the stapler and pry it open.",
            "risk": 10, # Highest risk: High chance of uncontrolled tear of the cecum.
            "effectiveness": 2 # Low effectiveness: Likely to fail or cause a worse problem.
        },
        'C': {
            "description": "Extend the port of the stapler into a longer incision then use a right angle to get between the jaws of the stapler and pry it open.",
            "risk": 8,  # Very high risk: Prying is still dangerous, even with better access.
            "effectiveness": 4 # Low effectiveness: Still relies on a dangerous maneuver.
        },
        'D': {
            "description": "Extend the port of the stapler into a longer incision then complete an open appendectomy via that incision.",
            "risk": 1,  # Lowest risk: Standard, safe "conversion to open" procedure. Provides maximum control.
            "effectiveness": 10 # Highest effectiveness: The definitive and safest solution.
        },
        'E': {
            "description": "Make a midline incision and use a right angle to get between the jaws of the stapler and pry it open.",
            "risk": 9,  # Very high risk: Combines a dangerous maneuver (prying) with an unnecessarily large and misplaced incision.
            "effectiveness": 3 # Low effectiveness: Poor surgical choice.
        },
        'F': {
            "description": "Make a midline incision and complete the appendectomy via an open approach.",
            "risk": 4,  # Moderate risk: Midline incision is excessive and more morbid than necessary for an appendectomy.
            "effectiveness": 9 # Effective, but the approach is suboptimal compared to extending the port.
        }
    }

    # Find the best option by prioritizing the lowest risk, then the highest effectiveness.
    # We sort by risk (ascending) and then by effectiveness (descending, hence the negative sign).
    best_option_key = sorted(options.keys(), key=lambda k: (options[k]['risk'], -options[k]['effectiveness']))[0]
    best_option_details = options[best_option_key]

    print("Evaluating surgical options based on Risk (1-10) and Effectiveness (1-10).")
    print("The primary principle is patient safety, so the lowest risk is prioritized.")
    print("-" * 30)

    print(f"The chosen option is '{best_option_key}': {best_option_details['description']}")
    print("\nRationale:")
    print("This option represents the safest and most standard 'bailout' procedure. When a laparoscopic procedure becomes unsafe due to an instrument malfunction, converting to an open procedure provides the direct visualization and control needed to resolve the issue without causing further harm. Extending an existing port site is the most efficient and least traumatic way to accomplish this.")

    # Fulfilling the request to output numbers in a final equation format
    risk_val = best_option_details['risk']
    effectiveness_val = best_option_details['effectiveness']
    print("\nFinal Calculation:")
    print(f"Risk = {risk_val}")
    print(f"Effectiveness = {effectiveness_val}")
    print(f"The final decision is option {best_option_key}.")

if __name__ == '__main__':
    solve_surgical_dilemma()