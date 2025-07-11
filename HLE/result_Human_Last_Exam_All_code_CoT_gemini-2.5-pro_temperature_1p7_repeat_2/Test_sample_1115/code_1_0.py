# This script determines the polarization of a companion photon in an entangled pair
# based on the principle of conservation of angular momentum.

# 1. Define the known physical quantities.
# The initial and final total angular momentum of the atom (1S0 state) is 0.
J_total = 0
# One photon is measured with right-handed polarization. Its angular momentum
# projection (helicity) is +1.
helicity_photon_1 = 1

# 2. State the conservation of angular momentum principle for the emitted photons.
# The total angular momentum of the two photons must sum to the system's total change, which is 0.
# The equation is: (Helicity of Photon 1) + (Helicity of Photon 2) = 0
print("The conservation of angular momentum requires the photons' helicities to sum to zero.")

# 3. Present the equation with the known values to be solved.
# Per the instructions, we show each number in the final equation.
# The equation is: 1 + helicity_photon_2 = 0
print(f"The equation with the measured value is:")
print(f"{helicity_photon_1} + helicity_photon_2 = {J_total}")

# 4. Solve for the helicity of the second photon.
helicity_photon_2 = J_total - helicity_photon_1

# 5. Interpret the result. A helicity of -1 corresponds to left-handed polarization.
if helicity_photon_2 == -1:
    polarization_photon_2 = "Left-handed"
elif helicity_photon_2 == 1:
    polarization_photon_2 = "Right-handed"
else:
    polarization_photon_2 = "Undetermined"

print(f"\nSolving for the second photon, its helicity must be {helicity_photon_2}.")
print(f"Conclusion: A helicity of {helicity_photon_2} means the companion photon's polarization is {polarization_photon_2}.")