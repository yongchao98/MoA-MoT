import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
captured_output = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = captured_output


# --- Case Data ---
# This data represents the patient's profile and key findings from the clinical vignette.
patient_age = 53
chest_pain_score = 10
smoking_packs_per_day = 2
smoking_duration_years = 20

# Clinical findings represented as boolean flags (True for present, False for absent)
has_heavy_alcohol_history = True
has_imaging_wall_thickening = True
# The endoscopy was crucial: no visible mucosal abnormalities.
has_endoscopic_ulcers_or_plaques = False

# --- Diagnostic Scoring Logic ---
# We will assign points to each possible diagnosis based on the evidence.

diagnoses = {
    "A. Streptococcal esophagitis": 0,
    "B. Esophageal adenocarcinoma": 0,
    "C. Esophageal squamous cell carcinoma": 0,
    "D. GERD": 0,
    "E. Herpes esophagitis": 0
}

# Score for Esophageal Squamous Cell Carcinoma (SCC)
# Heavy smoking and alcohol are major risk factors for SCC.
if smoking_packs_per_day >= 2 and smoking_duration_years >= 20:
    diagnoses["C. Esophageal squamous cell carcinoma"] += 5
if has_heavy_alcohol_history:
    diagnoses["C. Esophageal squamous cell carcinoma"] += 5
# Imaging showing wall thickening suggests an infiltrative process like a tumor.
if has_imaging_wall_thickening:
    diagnoses["C. Esophageal squamous cell carcinoma"] += 4
# A normal endoscopy with positive imaging is a classic sign of a submucosal tumor.
if not has_endoscopic_ulcers_or_plaques and has_imaging_wall_thickening:
    diagnoses["C. Esophageal squamous cell carcinoma"] += 5

# Score for other diagnoses (they are less likely)
# Infectious esophagitis (Strep/Herpes) would typically show ulcers or plaques.
if not has_endoscopic_ulcers_or_plaques:
    diagnoses["A. Streptococcal esophagitis"] += 1
    diagnoses["E. Herpes esophagitis"] += 1
# Adenocarcinoma is less strongly linked to these specific risk factors than SCC.
diagnoses["B. Esophageal adenocarcinoma"] += 2
# Severe GERD would likely show endoscopic changes, and it doesn't explain the wall thickening.
if not has_endoscopic_ulcers_or_plaques and not has_imaging_wall_thickening:
    diagnoses["D. GERD"] += 2
else:
    diagnoses["D. GERD"] += 0

# --- Print Final Analysis ---
# This section outputs the rationale and results.

# Fulfilling the requirement to print numbers from the prompt in an equation-like format
print("Diagnostic Analysis based on Patient Data:")
print(f"Patient Age ({patient_age}) + Smoking History ({smoking_packs_per_day} packs * {smoking_duration_years} years) + Pain Score ({chest_pain_score}/10) -> Diagnosis")
print("-" * 50)
print("Likelihood scores for each diagnosis:")
for diagnosis, score in diagnoses.items():
    print(f"{diagnosis}: {score}")

print("-" * 50)
most_likely = max(diagnoses, key=diagnoses.get)
print(f"Conclusion: The combination of major risk factors (smoking, alcohol) with imaging showing wall thickening\nand a 'normal' endoscopy is highly suggestive of a submucosally infiltrating tumor.")
print(f"The most likely diagnosis is: {most_likely}")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()

# Print the final result to the user
print(output_string)