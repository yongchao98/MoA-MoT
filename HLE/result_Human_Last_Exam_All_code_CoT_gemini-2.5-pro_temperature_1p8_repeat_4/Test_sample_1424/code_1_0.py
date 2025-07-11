import cmath

def calculate_bare_greens_function(omega, epsilon_k, eta):
  """
  Calculates the bare retarded Green's function G_0.

  The functional dependence of the bare Green's function G_0 on the
  single-particle energy eigenvalue epsilon_k is given by the formula:
  G_0(omega, epsilon_k) = 1 / (omega - epsilon_k + i*eta)
  where omega is the energy/frequency, and eta is a small positive infinitesimal.

  Args:
    omega (float): The energy or frequency variable.
    epsilon_k (float): The single-particle energy eigenvalue.
    eta (float): A small positive value to handle the pole.

  Returns:
    complex: The value of the bare Green's function.
  """
  # The denominator is a complex number
  denominator = (omega - epsilon_k) + 1j * eta
  # The Green's function is the inverse of the denominator
  g0 = 1 / denominator
  return g0

def main():
  """
  Main function to demonstrate the calculation.
  """
  # --- Example Values ---
  # Energy/frequency at which to evaluate the function
  omega = 2.5
  # Energy eigenvalue of a single-particle state
  epsilon_k = 2.0
  # Small positive infinitesimal
  eta = 0.1

  # Calculate the Green's function
  result = calculate_bare_greens_function(omega, epsilon_k, eta)

  # --- Output ---
  print("The functional dependence of the bare Green's function G_0 on the single-particle energy eigenvalue epsilon_k is:")
  print("G_0 = 1 / (omega - epsilon_k + i*eta)\n")

  print("For the example values:")
  print(f"  omega     = {omega}")
  print(f"  epsilon_k = {epsilon_k}")
  print(f"  eta       = {eta}\n")

  # Print the equation with the numbers substituted in, as requested.
  # We show each number explicitly in the final equation string.
  print("Final Equation:")
  print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
  
  # Print the final numerical result
  print(f"\nResult: G_0 = {result}")

if __name__ == "__main__":
  main()