import sys
import io

# Redirect stdout to capture print output for the final response
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Clinical Data ---
patient_age = 57
symptoms = ["dyspnea", "chronic cough", "acid reflux"]
history = "COPD"
imaging_finding = "mass of the vertebrae"
lab_finding_creatinine = 2.1

# --- Diagnostic Reasoning Model ---

# Define weights for key clinical findings. The vertebral mass is the most
# specific and high-impact finding, so it gets the highest weight.
weights = {
    'explains_vertebral_mass': 10,
    'explains_respiratory_symptoms': 3,
    'fits_patient_profile': 2
}

# Define how well each diagnosis explains the findings (1 for yes, 0 for no)
diagnoses_explanation = {
    'A. Aspiration pneumonitis': {'explains_vertebral_mass': 0, 'explains_respiratory_symptoms': 1, 'fits_patient_profile': 1},
    'B. Aspiration pneumonia':    {'explains_vertebral_mass': 0, 'explains_respiratory_symptoms': 1, 'fits_patient_profile': 1},
    'C. Achalasia':               {'explains_vertebral_mass': 0, 'explains_respiratory_symptoms': 0, 'fits_patient_profile': 0},
    'D. Adenocarcinoma':          {'explains_vertebral_mass': 1, 'explains_respiratory_symptoms': 1, 'fits_patient_profile': 1},
    'E. COPD':                    {'explains_vertebral_mass': 0, 'explains_respiratory_symptoms': 1, 'fits_patient_profile': 1}
}

# Calculate scores for each diagnosis
scores = {}
for diagnosis, explanations in diagnoses_explanation.items():
    score = (explanations['explains_vertebral_mass'] * weights['explains_vertebral_mass'] +
             explanations['explains_respiratory_symptoms'] * weights['explains_respiratory_symptoms'] +
             explanations['fits_patient_profile'] * weights['fits_patient_profile'])
    scores[diagnosis] = score

# Determine the most likely diagnosis
most_likely_diagnosis = max(scores, key=scores.get)
best_score_details = diagnoses_explanation[most_likely_diagnosis]
best_score_value = scores[most_likely_diagnosis]

# --- Output Generation ---
print("Analyzing the patient's case to find the most unifying diagnosis.")
print(f"Key Finding 1 (Imaging): {imaging_finding}")
print(f"Key Finding 2 (Symptoms): {', '.join(symptoms)}")
print(f"Key Finding 3 (History/Age): {patient_age}-year-old with {history}")
print("-" * 30)
print(f"The most likely diagnosis is: {most_likely_diagnosis}")
print("This conclusion is based on a scoring model that prioritizes the diagnosis's ability to explain the most critical findings.")
print("\nFinal Equation for the best diagnosis:")

# Print the equation with each number as requested
w_mass = weights['explains_vertebral_mass']
f_mass = best_score_details['explains_vertebral_mass']
w_resp = weights['explains_respiratory_symptoms']
f_resp = best_score_details['explains_respiratory_symptoms']
w_prof = weights['fits_patient_profile']
f_prof = best_score_details['fits_patient_profile']

print(f"Score = (Explains Mass? {f_mass} * Weight {w_mass}) + (Explains Resp. Symptoms? {f_resp} * Weight {w_resp}) + (Fits Profile? {f_prof} * Weight {w_prof})")
print(f"Score = ({f_mass * w_mass}) + ({f_resp * w_resp}) + ({f_prof * w_prof}) = {best_score_value}")
print("\nRationale: Adenocarcinoma (a type of cancer) is the only option that explains the vertebral mass, which is likely a metastasis from a primary tumor (e.g., in the lung) that is also causing the respiratory symptoms.")

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())