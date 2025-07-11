import math

def calculate_height(r, theta_deg):
  """
  Calculates the height of the cylinder based on its radius and the inscribed angle.

  Args:
    r: The radius of the cylinder.
    theta_deg: The inscribed angle in degrees.
  """
  # Convert angle from degrees to radians for math.tan()
  theta_rad = math.radians(theta_deg)

  # Calculate the height using the derived formula
  h = 2 * r * math.tan(theta_rad)

  # Print the formula with the given values
  print("The relationship is given by the formula: h = 2 * r * tan(theta)")
  print(f"Given r = {r} and theta = {theta_deg} degrees ({theta_rad:.4f} radians):")
  print(f"h = 2 * {r} * tan({theta_deg})")
  print(f"h = {h:.4f}")

# Example usage with r=10 and theta=30 degrees
example_r = 10
example_theta_deg = 30
calculate_height(example_r, example_theta_deg)