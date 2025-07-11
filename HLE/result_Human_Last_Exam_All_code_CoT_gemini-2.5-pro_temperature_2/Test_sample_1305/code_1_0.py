import sys
import io

# Helper function to capture print output
def execute_and_capture_output(func):
    """Executes a function and captures its print output."""
    old_stdout = sys.stdout
    sys.stdout = captured_output = io.StringIO()
    func()
    sys.stdout = old_stdout
    return captured_output.getvalue()

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    This is a text-based analysis presented within a Python script as requested.
    There are no numerical equations to solve. The output is a reasoned medical conclusion.
    """
    print("Step 1: Analyzing the patient's clinical presentation.")
    print("-----------------------------------------------------")
    print("Patient: 31-year-old woman.")
    print("Symptoms: Progressive shortness of breath, fatigue, lightheadedness, and a history of faint childhood cyanosis.")
    print("Physical Exam Murmur: A systolic ejection murmur at the left upper sternal border (LUSB), which increases with inspiration.")
    print("ECG Findings: Left axis deviation (LAD) and signs of right ventricular hypertrophy (RVH).\n")

    print("Step 2: Synthesizing the key findings.")
    print("---------------------------------------")
    print("The core challenge is to find a single condition that explains:")
    print("1. A systolic ejection murmur at the pulmonic area (LUSB).")
    print("2. A right-sided cardiac event (murmur increases with inspiration).")
    print("3. Evidence of right ventricular overload (RVH).")
    print("4. A paradoxical Left Axis Deviation (LAD), as most conditions causing RVH lead to Right Axis Deviation.\n")

    print("Step 3: Evaluating the answer choices against the synthesized findings.")
    print("--------------------------------------------------------------------")
    print("A. Ebstein anomaly: Incorrect. Causes a holosystolic regurgitant murmur, not a systolic ejection murmur.")
    print("B. Patent ductus arteriosus: Incorrect. Causes a continuous 'machine-like' murmur.")
    print("C. Mitral valve prolapse: Incorrect. A left-sided event with a different murmur type (mid-systolic click).")
    print("D. Atrial septal defect (ASD): Plausible. An ASD causes increased flow across the pulmonic valve, creating a systolic ejection murmur at the LUSB that increases with inspiration, and leads to RVH. While a standard ASD (ostium secundum) causes right axis deviation, a less common but well-known type, the ostium primum ASD, is classically associated with Left Axis Deviation. This diagnosis fits all findings perfectly.")
    print("E. Hypertrophic cardiomyopathy: Incorrect. A left-sided condition associated with Left Ventricular Hypertrophy, not RVH.")
    print("F. Tricuspid stenosis: Incorrect. Causes a diastolic rumble, not a systolic murmur.")
    print("G. Ventricular septal defect: Less likely. The typical VSD murmur is holosystolic and at the lower sternal border.")
    print("\nConclusion: The combination of a pulmonic flow murmur, RVH, and LAD is the classic triad for an ostium primum atrial septal defect. This is the most likely diagnosis.")

# Execute the analysis and print the result.
# The user instruction "output each number in the final equation" is not applicable
# as there is no equation. The step-by-step reasoning is provided instead.
analysis_text = execute_and_capture_output(solve_medical_case)
print(analysis_text)
print("<<<D>>>")