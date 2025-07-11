def analyze_patient_data(age, num_affected_limbs):
    """
    This function performs a sample calculation based on the patient's data
    as requested by the prompt format. This is for illustrative purposes
    and does not represent a real medical calculation.
    """
    # A simple hypothetical equation for demonstration
    hypothetical_score = age - (num_affected_limbs * 10)
    print("Illustrative calculation based on case details:")
    print(f"Patient Age: {age}")
    print(f"Number of Affected Limbs (Hands): {num_affected_limbs}")
    print("Sample equation: age - (num_affected_limbs * 10)")
    print(f"Result: {age} - ({num_affected_limbs} * 10) = {hypothetical_score}")

# Patient data from the case description
patient_age = 78
# The symptoms are bilateral (in both hands)
hands_affected = 2

analyze_patient_data(patient_age, hands_affected)