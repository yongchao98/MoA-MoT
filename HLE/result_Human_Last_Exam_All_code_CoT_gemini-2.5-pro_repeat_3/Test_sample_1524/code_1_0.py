import sys

# This script explains the fundamental limit on the chemical potential for bosons
# undergoing Bose-Einstein condensation.

# Step 1: The Bose-Einstein Distribution
# The average number of bosons (n_i) in a single-particle quantum state 'i'
# with energy ε_i is given by the Bose-Einstein distribution formula.
print("Step 1: The Bose-Einstein Distribution")
print("The governing equation is the Bose-Einstein distribution for the occupation number n_i of a state with energy ε_i:")
# We represent the equation as a string. Note that there are no numbers to substitute in this fundamental equation.
equation = "n_i = 1 / (exp((ε_i - μ) / (k_B * T)) - 1)"
print(equation)
print("where:")
print("  μ   is the chemical potential")
print("  T   is the absolute temperature")
print("  k_B is the Boltzmann constant")
print("  ε_i is the energy of the i-th state")
print("-" * 30)

# Step 2: The Physical Constraint
# A core principle of quantum statistics is that occupation numbers cannot be negative.
print("Step 2: The Physical Constraint")
print("A physical system cannot have a negative number of particles in any state. Therefore, for all states i, n_i must be non-negative (n_i >= 0).")
print("-" * 30)

# Step 3: Deriving the Mathematical Limit on μ
# From the constraint n_i >= 0, we can derive a limit on μ.
print("Step 3: Deriving the Mathematical Limit")
print("For n_i to be positive, the denominator of the distribution must also be positive:")
print("exp((ε_i - μ) / (k_B * T)) - 1 > 0")
print("Rearranging this gives:")
print("exp((ε_i - μ) / (k_B * T)) > 1")
print("Since the natural logarithm ln(x) is an increasing function for x > 0, we can take the log of both sides without changing the inequality:")
print("(ε_i - μ) / (k_B * T) > ln(1)")
# The number 1 is part of this derivation.
print("(ε_i - μ) / (k_B * T) > 0")
print("Since temperature T and Boltzmann's constant k_B are positive, this leads to the fundamental condition:")
print("ε_i - μ > 0, which means μ < ε_i")
print("-" * 30)

# Step 4: The Ultimate Limit from the Ground State
# This condition must hold for ALL states, so it is most restrictive for the state with the lowest energy.
print("Step 4: The Ultimate Limit from the Ground State")
print("The condition μ < ε_i must be true for every single-particle state.")
print("This implies that μ must be less than the lowest possible energy, which is the ground-state energy, ε_0.")
print("Thus, the ultimate limit is: μ ≤ ε_0")
print("The chemical potential can get arbitrarily close to ε_0 but can never exceed it.")
print("-" * 30)

# Step 5: Evaluating the Answer Choices
print("Step 5: Evaluating the Answer Choices")
print("The question asks for the fundamental limit on μ. We have found this limit to be the ground-state energy, ε_0.")
print("Now we check which answer choice corresponds to this limit.")
print("Choice C is: 'μ must be equal to the chemical potential of a non-interacting Bose gas at zero temperature.'")
print("At zero temperature (T=0), a non-interacting Bose gas has all its particles in the ground state. In this specific case, the chemical potential becomes exactly equal to the ground-state energy: μ(T=0) = ε_0.")
print("Therefore, Choice C provides a correct physical definition for the limiting value ε_0.")
print("-" * 30)

# Final Answer
sys.stdout.write("<<<C>>>\n")