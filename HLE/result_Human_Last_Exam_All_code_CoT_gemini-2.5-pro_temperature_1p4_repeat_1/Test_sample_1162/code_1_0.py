import operator

# Patient's key clinical findings
patient_findings = {
    "aniridia": True,
    "abdominal_mass": True,
    "hypertension": True,
    "developmental_delay": True
}

# Scoring system for differential diagnosis
# Points are assigned based on the strength of association
scoring_weights = {
    "aniridia": {"Nephroblastoma": 5, "Neuroblastoma": 0, "Germ cell tumor": 0, "Astrocytoma": 0, "Ewing sarcoma": 0},
    "abdominal_mass": {"Nephroblastoma": 2, "Neuroblastoma": 2, "Germ cell tumor": 1, "Astrocytoma": 0, "Ewing sarcoma": 1},
    "hypertension": {"Nephroblastoma": 2, "Neuroblastoma": 2, "Germ cell tumor": 0, "Astrocytoma": 0, "Ewing sarcoma": 0},
    "developmental_delay": {"Nephroblastoma": 2, "Neuroblastoma": 0, "Germ cell tumor": 0, "Astrocytoma": 0, "Ewing sarcoma": 0}
}

# Initialize scores
diagnosis_scores = {
    "Nephroblastoma": 0,
    "Neuroblastoma": 0,
    "Germ cell tumor": 0,
    "Astrocytoma": 0,
    "Ewing sarcoma": 0
}

# Calculate scores based on patient findings
equation_components = []
for finding, present in patient_findings.items():
    if present:
        for diagnosis, score in scoring_weights[finding].items():
            diagnosis_scores[diagnosis] += score
        # Add to the equation string for the top diagnosis later
        if scoring_weights[finding]["Nephroblastoma"] > 0:
            equation_components.append(str(scoring_weights[finding]["Nephroblastoma"]))

# Determine the most likely diagnosis
most_likely_diagnosis = max(diagnosis_scores.items(), key=operator.itemgetter(1))[0]

print("Diagnostic Scoring Process:")
print("--------------------------")
for diagnosis, score in diagnosis_scores.items():
    print(f"Final Score for {diagnosis}: {score}")
print("--------------------------\n")


print(f"The clinical picture strongly points towards {most_likely_diagnosis}.")
print("This is due to the classic WAGR syndrome association (Wilms Tumor, Aniridia, Genitourinary anomalies, Retardation).")
print("The aniridia is a particularly specific finding.\n")

# Fulfilling the "equation" requirement for the most likely diagnosis
equation_str = " + ".join(equation_components)
final_score = diagnosis_scores[most_likely_diagnosis]
print(f"Scoring equation for {most_likely_diagnosis}: {equation_str} = {final_score}")
print("\nFinal Answer corresponds to Nephroblastoma (also known as Wilms tumor).")
