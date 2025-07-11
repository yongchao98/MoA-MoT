import math

# This script calculates the Bell parameter 'S' for the CHSH inequality.
# According to local realism (classical physics, or a dream simulation), |S| <= 2.
# According to quantum mechanics, |S| can be up to 2 * sqrt(2) ~= 2.828.
# A result > 2 would be strong evidence of being in an objective, quantum reality.

# We define the optimal angles for the polarizers (in radians)
# to achieve the maximum violation of the inequality.
angle_a = 0
angle_a_prime = math.pi / 4  # 45 degrees
angle_b = math.pi / 8  # 22.5 degrees
angle_b_prime = 3 * math.pi / 8 # 67.5 degrees

# The quantum mechanical expectation value for the correlation between two polarizers
# with an angle 'theta' between them is given by -cos(2 * theta).
# E(a, b) is the correlation for detectors set at angles 'a' and 'b'.
def expectation_value(angle1, angle2):
    """Calculates the correlation E(a,b) = -cos(2 * (b - a))."""
    return -math.cos(2 * (angle2 - angle1))

# Calculate the four correlation terms needed for the S value.
E_ab = expectation_value(angle_a, angle_b)
E_ab_prime = expectation_value(angle_a, angle_b_prime)
E_a_prime_b = expectation_value(angle_a_prime, angle_b)
E_a_prime_b_prime = expectation_value(angle_a_prime, angle_b_prime)

# The CHSH parameter 'S' is calculated using the formula:
# S = E(a, b) - E(a, b') + E(a', b) + E(a', b')
S = E_ab - E_ab_prime + E_a_prime_b + E_a_prime_b_prime

print("Calculating the CHSH parameter 'S' using quantum mechanics:")
print(f"S = E(a,b) - E(a,b') + E(a',b) + E(a',b')")
print(f"S = ({E_ab:.4f}) - ({E_ab_prime:.4f}) + ({E_a_prime_b:.4f}) + ({E_a_prime_b_prime:.4f})")
print(f"Final calculated S value: {S:.4f}")
print("\n---")
print("Comparison:")
print(f"Classical/Dream Limit: |S| <= 2")
print(f"Quantum Mechanical Result: |S| = {abs(S):.4f}")

if abs(S) > 2:
    print("Conclusion: The result violates the classical limit, indicating an objective, quantum reality.")
else:
    print("Conclusion: The result is consistent with a classical simulation (a dream).")
