def solve_medical_question():
    """
    This script presents the solution to the medical question based on the case study.
    """
    # Patient's relevant numbers from the case study
    systolic_bp = 130
    diastolic_bp = 70

    # The correct mechanism
    mechanism = "Through direct vasodilation"
    explanation = "Magnesium acts as a natural calcium antagonist (a calcium channel blocker), which helps relax vascular smooth muscle, leading to vasodilation and a decrease in blood pressure."

    # Outputting the 'equation' as a statement connecting the patient's numbers to the answer
    print(f"The patient's blood pressure is {systolic_bp}/{diastolic_bp} mm Hg.")
    print(f"Magnesium supplementation can help lower blood pressure via the following primary mechanism:")
    print(f"Answer: {mechanism}")
    print(f"Explanation: {explanation}")

solve_medical_question()