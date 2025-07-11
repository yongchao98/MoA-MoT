import pandas as pd

def solve_surgical_dilemma():
    """
    Analyzes the surgical options for a stuck laparoscopic stapler during an appendectomy
    and determines the best course of action.
    """
    # Define the options and their scores based on key surgical principles.
    # Scores are on a scale of 1-10.
    # Safety: Higher is better.
    # Effectiveness: Higher is better.
    # Invasiveness: Lower is better.
    data = {
        "A": {
            "Description": "Staple and resect part of the cecum.",
            "Safety": 1, "Effectiveness": 2, "Invasiveness": 7
        },
        "B": {
            "Description": "Pry open laparoscopically.",
            "Safety": 2, "Effectiveness": 3, "Invasiveness": 3
        },
        "C": {
            "Description": "Extend port and pry open.",
            "Safety": 4, "Effectiveness": 5, "Invasiveness": 6
        },
        "D": {
            "Description": "Extend port and convert to open appendectomy.",
            "Safety": 9, "Effectiveness": 10, "Invasiveness": 6
        },
        "E": {
            "Description": "Make midline incision and pry open.",
            "Safety": 4, "Effectiveness": 5, "Invasiveness": 10
        },
        "F": {
            "Description": "Make midline incision and convert to open appendectomy.",
            "Safety": 8, "Effectiveness": 10, "Invasiveness": 10
        }
    }

    # Convert to a pandas DataFrame for easier manipulation
    df = pd.DataFrame.from_dict(data, orient='index')

    # Calculate a weighted score. Safety is the most important factor.
    # Score = (Safety * 2) + (Effectiveness * 1.5) - (Invasiveness * 1)
    df['Total_Score'] = (df['Safety'] * 2) + (df['Effectiveness'] * 1.5) - df['Invasiveness']

    # Find the best option
    best_option_key = df['Total_Score'].idxmax()
    best_option_details = df.loc[best_option_key]

    print("Analysis of Surgical Options:\n")
    print(df)
    print("\n-------------------------------------------------\n")
    print("Conclusion:\n")
    print("The primary goal is patient safety. The laparoscopic approach has failed due to an equipment malfunction.")
    print("Continuing with risky laparoscopic maneuvers (A, B) or making an unnecessarily large incision (E, F) is not ideal.")
    print("The safest and most effective approach is to convert the procedure to an open one through a minimal but adequate incision.")
    print("This allows for direct visualization, control of the anatomy, and safe completion of the appendectomy.")
    print("\nBased on the analysis, the best option is:")
    print(f"Option {best_option_key}: {best_option_details['Description']}")
    print(f"This option has the highest calculated score of {best_option_details['Total_Score']:.1f}, balancing safety, effectiveness, and appropriate level of invasiveness.")

solve_surgical_dilemma()