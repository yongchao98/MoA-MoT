import math

def calculate_completeness_bound(b, d):
  """
  Calculates the theoretical completeness bound based on the formula ln(b)/d.

  Args:
    b (int): The branching factor of the source ontology O₁.
    d (int): The depth of the source ontology O₁.
  """
  if b <= 0 or d <= 0:
    print("Branching factor and depth must be positive numbers.")
    return

  # Calculate the theoretical completeness bound
  completeness_bound = math.log(b) / d

  # Print the explanation and the result, including the numbers from the equation
  print("This scenario describes a loss of information at each level of the hierarchy due to the conflict between 'mergeable siblings' and a 'strict categorical structure'.")
  print("The theoretical completeness of the mapping is bounded by a measure related to the information preserved across the hierarchy.")
  print(f"The formula for this bound is ln(b)/d.")
  print(f"Given b = {b} and d = {d}, the calculation is:")
  print(f"ln({b}) / {d} = {completeness_bound:.4f}")
  print(f"\nThis value, approximately {completeness_bound:.2%}, represents the upper limit on the proportion of structural information that can be retained in the mapping.")

# Parameters from the problem description
branching_factor = 3
depth = 4

# Execute the calculation and print the results
calculate_completeness_bound(branching_factor, depth)