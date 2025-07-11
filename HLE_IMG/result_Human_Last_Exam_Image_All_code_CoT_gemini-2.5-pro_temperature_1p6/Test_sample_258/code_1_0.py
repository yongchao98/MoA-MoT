import math

def calculate_cp_length(ab_length):
  """
  Calculates the length of CP based on the length of AB.
  
  The derivation shows that CP = AB * φ, where φ is the golden ratio.
  φ = (1 + sqrt(5)) / 2
  
  Args:
    ab_length: The length of the segment AB.
    
  Returns:
    The length of the segment CP.
  """
  
  # The golden ratio
  phi = (1 + math.sqrt(5)) / 2
  
  # Calculate the length of CP
  cp_length = ab_length * phi
  
  # Print the equation with the given values
  print(f"Given AB = {ab_length}")
  print("The relationship between CP and AB is: CP = AB * ((sqrt(5) + 1) / 2)")
  print(f"So, CP = {ab_length} * (({math.sqrt(5):.4f} + 1) / 2)")
  print(f"CP = {ab_length} * ({(math.sqrt(5) + 1):.4f} / 2)")
  print(f"CP = {ab_length} * {((math.sqrt(5) + 1) / 2):.4f}")
  print(f"The calculated length of CP is: {cp_length:.4f}")

# The given length of AB
ab_length = 18468
calculate_cp_length(ab_length)

# The final answer in the required format
final_answer = ab_length * (1 + math.sqrt(5)) / 2
# print(f"\n<<<{final_answer:.4f}>>>")