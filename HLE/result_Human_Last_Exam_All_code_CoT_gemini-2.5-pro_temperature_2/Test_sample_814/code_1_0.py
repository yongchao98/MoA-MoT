import pandas as pd

def solve_clinical_case():
    """
    Analyzes a clinical case of suspected fibromyalgia and determines the best treatment option.
    """
    # Patient Profile
    patient_symptoms = {
        'Pain': 'Widespread, >1 year',
        'Energy': 'Extreme fatigue',
        'Mood': 'Anxiety and Depression',
        'Sleep': 'Sleep issues, Restless Leg Syndrome',
        'Cognition': 'Diminished cognitive ability',
        'Neurological': 'Paraesthesia',
        'Current Meds': 'Ibuprofen (minimal relief)',
        'Ruled Out': 'Thyroid disease, Rheumatoid Arthritis, Lupus'
    }

    # Answer Choices
    choices = {
        'A': 'Duloxetine+Gabapentin',
        'B': 'Gabapentin',
        'C': 'Duloxetine',
        'D': 'cyclobenzaprine',
        'E': 'Duloxetine+acetaminophen',
        'F': 'Duloxetine+cyclobenzaprine'
    }

    print("Step 1: Summarize Patient Presentation and Likely Diagnosis")
    print("---------------------------------------------------------")
    print("The patient exhibits classic signs of Fibromyalgia, including widespread pain, fatigue, sleep/mood issues, and cognitive difficulties, with other major causes ruled out.\n")

    print("Step 2: Evaluate How Each Treatment Option Addresses the Patient's Symptoms")
    print("--------------------------------------------------------------------------")
    print("We will evaluate each option based on its effectiveness for Pain, Mood, Sleep, and Restless Legs/Paraesthesia.\n")

    # Suitability scores (out of 10) for each key symptom area
    evaluation = {
        'Pain':          {'A': 9, 'B': 7, 'C': 7, 'D': 4, 'E': 7, 'F': 8},
        'Mood':          {'A': 8, 'B': 2, 'C': 8, 'D': 1, 'E': 8, 'F': 8},
        'Sleep':         {'A': 9, 'B': 8, 'C': 5, 'D': 7, 'E': 5, 'F': 8},
        'RLS/Paraesth.': {'A': 9, 'B': 9, 'C': 3, 'D': 1, 'E': 3, 'F': 3}
    }
    df = pd.DataFrame(evaluation)
    df['Total_Score'] = df.sum(axis=1)

    print("Evaluation Matrix (Score out of 10 for each category):")
    print(df)
    print("\n")

    print("Step 3: Analyze the Total Scores and Rationale")
    print("---------------------------------------------")
    best_choice = df['Total_Score'].idxmax()
    
    print(f"Choice A (Duloxetine+Gabapentin) Total Score = {df.loc['A', 'Pain']} + {df.loc['A', 'Mood']} + {df.loc['A', 'Sleep']} + {df.loc['A', 'RLS/Paraesth.']} = {int(df.loc['A', 'Total_Score'])}")
    print("   - Rationale: Excellent synergy. Duloxetine targets pain and mood. Gabapentin adds potent pain relief, improves sleep, and is a first-line treatment for restless legs and paresthesia. This combination addresses the full spectrum of the patient's severe symptoms.")

    print(f"Choice B (Gabapentin) Total Score = {df.loc['B', 'Pain']} + {df.loc['B', 'Mood']} + {df.loc['B', 'Sleep']} + {df.loc['B', 'RLS/Paraesth.']} = {int(df.loc['B', 'Total_Score'])}")
    print("   - Rationale: Good for pain, sleep, and RLS, but does not adequately address the significant anxiety and depression.")

    print(f"Choice C (Duloxetine) Total Score = {df.loc['C', 'Pain']} + {df.loc['C', 'Mood']} + {df.loc['C', 'Sleep']} + {df.loc['C', 'RLS/Paraesth.']} = {int(df.loc['C', 'Total_Score'])}")
    print("   - Rationale: An excellent first-line monotherapy for pain and mood, but may not be sufficient for the severe sleep issues and RLS.")
    
    print(f"Choice D (cyclobenzaprine) Total Score = {df.loc['D', 'Pain']} + {df.loc['D', 'Mood']} + {df.loc['D', 'Sleep']} + {df.loc['D', 'RLS/Paraesth.']} = {int(df.loc['D', 'Total_Score'])}")
    print("   - Rationale: Mainly helps with sleep. Insufficient for comprehensive treatment.")
    
    print(f"Choice E (Duloxetine+acetaminophen) Total Score = {df.loc['E', 'Pain']} + {df.loc['E', 'Mood']} + {df.loc['E', 'Sleep']} + {df.loc['E', 'RLS/Paraesth.']} = {int(df.loc['E', 'Total_Score'])}")
    print("   - Rationale: Adding acetaminophen offers little to no benefit for this type of pain.")

    print(f"Choice F (Duloxetine+cyclobenzaprine) Total Score = {df.loc['F', 'Pain']} + {df.loc['F', 'Mood']} + {df.loc['F', 'Sleep']} + {df.loc['F', 'RLS/Paraesth.']} = {int(df.loc['F', 'Total_Score'])}")
    print("   - Rationale: A decent option for pain, mood, and sleep, but Gabapentin is superior for the RLS/paresthesia component.")

    print("\nStep 4: Final Conclusion")
    print("--------------------------")
    print(f"Based on the analysis, the combination in option '{best_choice}' provides the most comprehensive coverage for the patient's symptoms.")

solve_clinical_case()