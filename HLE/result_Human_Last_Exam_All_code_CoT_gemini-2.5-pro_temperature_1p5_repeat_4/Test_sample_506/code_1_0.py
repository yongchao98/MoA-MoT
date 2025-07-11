import sys
import io

# Buffer to capture output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Medication Recommendations ---
# This patient has resistant hypertension, diabetes, and is African American.
# Guidelines suggest a combination of an ACEi/ARB, a CCB, and a diuretic.
# The patient has many diuretics and other specific drugs she cannot take.

# 1. ACE Inhibitor: Essential for a patient with diabetes for renal protection.
medication_1 = "Lisinopril"

# 2. Calcium Channel Blocker: Excellent choice for African American patients.
# The patient cannot take verapamil (non-DHP), so a DHP-CCB is appropriate.
medication_2 = "Amlodipine"

# 3. Add-on for Resistant Hypertension: Guidelines recommend an MRA like spironolactone.
# Most standard diuretics are on the exclusion list, but spironolactone is not.
medication_3 = "Spironolactone"

# --- Output the recommendation ---
print("Based on the patient's profile, the following three medications are recommended to maximize hypertension treatment:")
print(f"1. {medication_1}")
print(f"2. {medication_2}")
print(f"3. {medication_3}")
print("\nRationale: This regimen combines an ACE inhibitor (essential for diabetes), a calcium channel blocker (highly effective in African American patients), and an aldosterone antagonist (specifically recommended for resistant hypertension when other diuretics are not an option).")

# --- Final processing to print captured output ---
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)

# The final answer format required by the prompt
final_answer = f"<<<{medication_1}, {medication_2}, {medication_3}>>>"
# This final_answer line is for the system and would not be printed in a real script.
# print(final_answer) # This would be uncommented if the execution environment needed to see it.