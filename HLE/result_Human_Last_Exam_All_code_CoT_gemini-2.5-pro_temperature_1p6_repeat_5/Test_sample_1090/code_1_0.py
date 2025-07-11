def calculate_knee_contracture():
    """
    Calculates the knee flexion contracture based on the patient's examination findings.
    """
    full_extension_degrees = 180
    max_extension_achieved = 165
    
    contracture_degrees = full_extension_degrees - max_extension_achieved
    
    print(f"The patient's inability to walk is partially due to a mechanical block in his left knee.")
    print(f"A fully extended knee is {full_extension_degrees} degrees.")
    print(f"The patient's knee can only extend to {max_extension_achieved} degrees.")
    print(f"Calculation: {full_extension_degrees} - {max_extension_achieved} = {contracture_degrees}")
    print(f"Result: The patient has a {contracture_degrees}-degree flexion contracture of the left knee.")

calculate_knee_contracture()