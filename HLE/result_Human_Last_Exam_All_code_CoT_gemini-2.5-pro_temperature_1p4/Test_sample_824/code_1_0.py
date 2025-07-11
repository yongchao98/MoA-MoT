def calculate_complement_drop(baseline_c3, percentage_drop):
    """
    Calculates the new complement level after a specified percentage drop.

    Args:
        baseline_c3 (float): The initial, normal C3 level.
        percentage_drop (float): The percentage by which the C3 level dropped.
    """
    drop_amount = baseline_c3 * (percentage_drop / 100)
    final_c3 = baseline_c3 - drop_amount

    print(f"Patient's Baseline C3 Level: {baseline_c3} mg/dL (Normal range: 90-180 mg/dL)")
    print(f"Percentage Drop during Renal Flare: {percentage_drop}%")
    print(f"Amount of C3 consumed: {baseline_c3} * ({percentage_drop}/100) = {drop_amount:.2f} mg/dL")
    print(f"Patient's C3 Level during Acute Flare: {baseline_c3} - {drop_amount:.2f} = {final_c3:.2f} mg/dL")
    print("\nThis significant drop in C3, a key complement protein, indicates high consumption due to immune complex deposition in the kidneys, which is the direct cause of severe lupus nephritis.")

# --- Patient's Hypothetical Scenario ---
# A typical normal baseline C3 level.
patient_baseline_c3 = 120.0
# A severe drop indicating an acute flare.
patient_percentage_drop = 65.0

calculate_complement_drop(patient_baseline_c3, patient_percentage_drop)