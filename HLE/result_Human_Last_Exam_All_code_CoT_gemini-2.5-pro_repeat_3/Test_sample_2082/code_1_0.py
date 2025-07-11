import sys

def solve_medical_case():
    # Patient data from the case study
    weight_kg = 65.77
    height_m = 1.68

    # --- Step 1: Perform a relevant calculation (BMI) ---
    # The equation is BMI = weight (kg) / (height (m))^2
    bmi = weight_kg / (height_m ** 2)

    print("Part 1: Clinical Calculation from Patient Data")
    print("---------------------------------------------")
    print("Calculating the patient's Body Mass Index (BMI).")
    print(f"The patient's weight is {weight_kg} kg and height is {height_m} m.")
    # Here we output each number in the final equation as requested
    print(f"BMI = {weight_kg} / ({height_m} * {height_m})")
    print(f"The patient's BMI is {bmi:.2f}, which falls within the normal weight range.\n")

    # --- Step 2: Explain the answer to the multiple-choice question ---
    print("Part 2: Answering the Question about Magnesium")
    print("---------------------------------------------")
    print("Question: By which mechanism can magnesium supplementation help lower blood pressure?")
    print("\nExplanation:")
    print("Magnesium helps lower blood pressure through several mechanisms, but a primary one is its role as a natural calcium channel blocker.")
    print("1. Calcium ions are responsible for the contraction of smooth muscle cells in the walls of blood vessels. When these muscles contract, the vessels narrow (vasoconstriction), increasing blood pressure.")
    print("2. Magnesium directly competes with calcium for entry into these cells. By blocking calcium's entry, magnesium promotes the relaxation of the vascular smooth muscle.")
    print("3. This relaxation leads to the widening of blood vessels (vasodilation), which reduces resistance to blood flow and thereby lowers systemic blood pressure.")
    print("Therefore, magnesium helps lower blood pressure through direct vasodilation.")

# Execute the function
solve_medical_case()