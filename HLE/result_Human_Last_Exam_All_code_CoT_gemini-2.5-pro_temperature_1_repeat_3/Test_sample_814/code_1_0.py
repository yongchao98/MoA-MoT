import sys
import io

# Block stdout to avoid printing inside the logic, will print everything at the end.
original_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The patient presents with a classic picture of Fibromyalgia.
# We need to choose a treatment that covers the wide range of symptoms.
# Let's model this decision with a scoring system.

# 1. List patient's symptoms
symptoms = {
    "widespread_pain": "Widespread Pain",
    "fatigue": "Fatigue",
    "anxiety": "Anxiety",
    "depression": "Depression",
    "sleep_issues": "Sleep Issues",
    "cognitive_diminution": "Cognitive Diminution ('Fibro Fog')",
    "restless_leg": "Restless Leg Syndrome",
    "paraesthesia": "Paraesthesia (Tingling)"
}

# 2. Assign efficacy scores to each drug for each symptom (0=None, 1=Low, 2=Moderate, 3=High)
drug_efficacy = {
    "Duloxetine": {
        "widespread_pain": 3, "fatigue": 2, "anxiety": 3, "depression": 3,
        "sleep_issues": 1, "cognitive_diminution": 2, "restless_leg": 0, "paraesthesia": 1
    },
    "Gabapentin": {
        "widespread_pain": 2, "fatigue": 0, "anxiety": 1, "depression": 0,
        "sleep_issues": 3, "cognitive_diminution": 0, "restless_leg": 3, "paraesthesia": 3
    },
    "cyclobenzaprine": {
        "widespread_pain": 2, "fatigue": -1, "anxiety": 0, "depression": 0,
        "sleep_issues": 3, "cognitive_diminution": -1, "restless_leg": 0, "paraesthesia": 0
    },
    "acetaminophen": {
        "widespread_pain": 1, "fatigue": 0, "anxiety": 0, "depression": 0,
        "sleep_issues": 0, "cognitive_diminution": 0, "restless_leg": 0, "paraesthesia": 0
    }
}

# 3. Define the answer choices as combinations of drugs
answer_choices = {
    "A": ["Duloxetine", "Gabapentin"],
    "B": ["Gabapentin"],
    "C": ["Duloxetine"],
    "D": ["cyclobenzaprine"],
    "E": ["Duloxetine", "acetaminophen"],
    "F": ["Duloxetine", "cyclobenzaprine"]
}

# 4. Calculate the total efficacy score for each choice
results = {}
calculation_steps = []

for choice, drugs in answer_choices.items():
    total_score = 0
    equation = []
    
    choice_str = f"Option {choice} ({' + '.join(drugs)})"
    calculation_steps.append(f"\n--- Evaluating {choice_str} ---")
    
    for key, name in symptoms.items():
        # For a combination, the efficacy for a symptom is the max of the individual drug efficacies
        symptom_score = max(drug_efficacy.get(drug, {}).get(key, 0) for drug in drugs)
        total_score += symptom_score
        equation.append(str(symptom_score))
    
    results[choice] = total_score
    final_equation = " + ".join(equation)
    calculation_steps.append(f"Calculation for {choice}: {final_equation} = {total_score}")

# 5. Find the best option
best_choice = max(results, key=results.get)

# Restore stdout and print the captured output
sys.stdout = original_stdout
print("Analysis of Treatment Options for Fibromyalgia Symptoms:\n")
print("\n".join(calculation_steps))

print("\n--- Summary of Scores ---")
for choice, score in sorted(results.items(), key=lambda item: item[1], reverse=True):
    print(f"Option {choice}: Total Score = {score}")

print(f"\nConclusion: Option {best_choice} provides the most comprehensive symptom coverage.")
print("This combination effectively targets the core fibromyalgia symptoms of pain and mood (Duloxetine) while also addressing the common co-occurring neuropathic symptoms like paraesthesia and restless leg syndrome, as well as sleep issues (Gabapentin).")

<<<A>>>