import sys
import io

# Capture original stdout to restore it later if needed, although not strictly necessary here.
original_stdout = sys.stdout
# Create a string buffer to hold the output.
output_buffer = io.StringIO()
# Redirect stdout to the buffer.
sys.stdout = output_buffer

# --- Diagnostic Logic Starts Here ---

# Step 1: Define the patient's clinical profile
patient_age = 18
patient_sex = "female"
height_percentile_lower = 10
height_percentile_upper = 15

presenting_symptoms = {
    "reproductive": ["amenorrhea", "infertility"],
    "growth": ["short stature"],
    "systemic": ["fatigue and shortness of breath on exertion", "occasional elevated blood pressure"]
}

exam_and_test_findings = {
    "ultrasound": "ovarian dysgenesis, normal kidneys",
    "genetics": "normal chromosomal complement (karyotype)"
}

print("--- Patient Data Analysis ---")
print(f"Patient is a {patient_age}-year-old {patient_sex}.")
print(f"Key Findings:")
print(f"  - Growth: Persistent short stature (in the {height_percentile_lower}th - {height_percentile_upper}th percentile).")
print(f"  - Ovarian function: Ovarian dysgenesis leading to amenorrhea and infertility.")
print(f"  - Other signs: Symptoms suggestive of cardiac strain and occasional hypertension.")
print(f"  - Karyotype: Reported as normal, which for a phenotypic female implies 46,XX.")
print("-" * 30 + "\n")


print("--- Diagnostic Reasoning ---")

# Step 2: Main diagnostic logic based on the combination of findings
is_turner_phenotype = ("ovarian dysgenesis" in exam_and_test_findings["ultrasound"] and
                         "short stature" in presenting_symptoms["growth"])

is_classic_turner_syndrome = (exam_and_test_findings["genetics"] != "normal chromosomal complement (karyotype)")

if is_turner_phenotype:
    print("Logic: The patient presents with the classic features of Turner Syndrome (short stature and ovarian dysgenesis).")
    if not is_classic_turner_syndrome:
        print("Logic: However, a diagnosis of classic Turner Syndrome (45,X0) is excluded by the normal karyotype result.")
        print("Logic: This situation points to a molecular abnormality that causes a 'Turner-like' phenotype without being a full chromosome loss.")
        print("\nConclusion:")
        print("The features of Turner Syndrome result from having only one functional copy of key genes on the X chromosome (haploinsufficiency).")
        print("  - The 'SHOX' gene on the short arm of the X chromosome is critical for height.")
        print("  - Other genes on the X chromosome are essential for ovarian development.")
        print("\nA single molecular event that can cause both short stature and ovarian failure in a 46,XX female is a deletion of genetic material from the short arm of one of the X chromosomes (Xp).")
        print("This 'Xp deletion' can be too small to be detected by standard karyotyping but is sufficient to cause the full clinical picture.")
        print("\nTherefore, the most likely molecular abnormality is:")
        final_answer = "Deletion of the short arm of the X chromosome (Xp deletion)"
        print(f">>> {final_answer}")
    else:
        final_answer = "Classic Turner Syndrome (45,X0)"
        print(f"This would be classic Turner Syndrome, but it contradicts the patient's test results.")
else:
    final_answer = "Diagnosis requires further investigation into other causes."
    print("The clinical picture is not a classic match for this specific line of reasoning.")

# --- Diagnostic Logic Ends Here ---

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
final_output = output_buffer.getvalue()

# Print the captured output
print(final_output.strip())