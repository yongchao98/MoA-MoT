import numpy as np

def calculate_l(b, n=101):
  """
  Calculates the value of l(b) based on the simplified formula Tr(K⁻¹).
  l(b) = (n + (n-2)*b^2) / (1-b^2)
  """
  b_squared = b**2
  numerator = n + (n - 2) * b_squared
  denominator = 1 - b_squared
  return numerator / denominator

def main():
  """
  Main function to compute the final expression.
  """
  # The value of b for the first term
  b1 = 0.5
  # The value of b for the second term
  b2 = -0.5

  # Calculate l(1/2)
  l_half = calculate_l(b1)
  
  # Calculate l(-1/2)
  l_minus_half = calculate_l(b2)

  # The expression to compute is 6 * (l(1/2) + l(-1/2))
  # Since l(b) is an even function, l(1/2) = l(-1/2)
  # So the expression is 12 * l(1/2)
  
  result = 6 * (l_half + l_minus_half)

  print("Step-by-step calculation:")
  print(f"n = 101")
  print(f"b₁ = {b1}, b₂ = {b2}")
  print(f"l(b) = (n + (n-2)*b²) / (1-b²)")
  
  # Showing calculation for l(1/2)
  num_half = 101 + 99 * (b1**2)
  den_half = 1 - (b1**2)
  print(f"l({b1}) = (101 + 99*({b1})²) / (1 - ({b1})²) = {num_half} / {den_half} = {l_half}")

  # Showing calculation for l(-1/2)
  num_mhalf = 101 + 99 * (b2**2)
  den_mhalf = 1 - (b2**2)
  print(f"l({b2}) = (101 + 99*({b2})²) / (1 - ({b2})²) = {num_mhalf} / {den_mhalf} = {l_minus_half}")
  
  print("\nFinal equation:")
  print(f"6 * (l({b1}) + l({b2})) = 6 * ({l_half} + {l_minus_half}) = {result}")
  
  # Final Answer
  print(f"\nThe computed value is: {result}")


if __name__ == "__main__":
  main()
