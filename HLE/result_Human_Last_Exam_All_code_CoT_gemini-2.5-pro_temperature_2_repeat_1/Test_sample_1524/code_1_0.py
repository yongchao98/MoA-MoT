import math

print("This script explains the fundamental limit on the chemical potential (μ) in Bose-Einstein condensation.")
print("-" * 70)

print("\nStep 1: The Bose-Einstein Distribution Function")
print("The average number of bosons, n(ε), in a state with energy ε is given by:")
print("n(ε) = 1 / (exp((ε - μ) / (k_B * T)) - 1)")
print("where μ is the chemical potential and ε_0 is the ground state energy (the lowest possible energy).\n")

print("Step 2: The Core Physical Constraint")
print("For the occupation number n(ε) to be positive, the denominator must be positive.")
print("This means: exp((ε - μ) / (k_B * T)) - 1 > 0")
print("Which simplifies to: exp((ε - μ) / (k_B * T)) > 1")
print("This inequality only holds if the exponent is positive: ε - μ > 0, or μ < ε.\n")

print("Step 3: The Universal Limit for μ")
print("The condition μ < ε must be true for ALL energy states ε.")
print("Therefore, μ must be less than the lowest possible energy, the ground state energy ε_0.")
print("The strict condition is: μ < ε_0.\n")

print("Step 4: Reaching the Limit for Condensation")
print("Bose-Einstein condensation occurs when a large number of particles occupy the ground state (ε_0).")
print("For n(ε_0) to be large, the denominator 'exp((ε_0 - μ) / (k_B*T)) - 1' must be very small.")
print("This happens as μ gets very close to ε_0 from below.")
print("In the condensed phase (at or below the critical temperature), the chemical potential becomes equal to the ground state energy in the thermodynamic limit: μ = ε_0.\n")

print("Step 5: Connecting to the Correct Answer Choice")
print("The chemical potential of a non-interacting Bose gas at absolute zero (T=0) is exactly the ground state energy, ε_0, because all particles are in that state.")
print("Therefore, the 'fundamental limit' that the chemical potential reaches for condensation is precisely the ground state energy ε_0.")
print("-" * 70)
print("The final equation describing the chemical potential in the condensed phase is:")
# The following lines demonstrate printing each "number" or symbol in the final equation.
chemical_potential = "μ"
equals_sign = "="
ground_state_energy = "ε_0"
print(f"Term 1 (Chemical Potential): {chemical_potential}")
print(f"Term 2 (Equality): {equals_sign}")
print(f"Term 3 (Ground State Energy): {ground_state_energy}")
print(f"Full Equation: {chemical_potential} {equals_sign} {ground_state_energy}")