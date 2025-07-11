import math

# Define arm segment lengths in cm
L1 = 40  # Shoulder to Elbow
L2 = 28  # Elbow to Wrist
L3 = 15  # Wrist to Hand
L4 = 10  # Hand to Finger

# The tightest collision constraint determines sin(alpha)
# dist(L1, L4) >= 1 cm gives (L2-L4)*sin(alpha) >= 1
# 18 * sin(alpha) >= 1  => sin(alpha) = 1/18
sin_alpha = 1.0 / (L2 - L4)

# Calculate cos(alpha) from sin(alpha)
# cos(alpha) = sqrt(1 - sin(alpha)^2)
cos_alpha = math.sqrt(1 - sin_alpha**2)

# The distance D of the fingertip from the shoulder in the chosen configuration is:
# D^2 = (L1-L3)^2 + (L2-L4)^2 - 2*(L1-L3)*(L2-L4)*cos(alpha)
term1 = (L1 - L3)**2
term2 = (L2 - L4)**2
term3 = 2 * (L1 - L3) * (L2 - L4) * cos_alpha

distance_sq = term1 + term2 - term3
distance = math.sqrt(distance_sq)

print("Step-by-step calculation of the final distance:")
print(f"1. Lengths: L1={L1}, L2={L2}, L3={L3}, L4={L4}")
print(f"2. The governing constraint is dist(L1, L4) >= 1 cm, which sets sin(alpha) = 1 / (L2 - L4) = 1/{int(L2-L4)}.")
print(f"3. sin(alpha) = {sin_alpha:.4f}")
print(f"4. cos(alpha) = sqrt(1 - sin(alpha)^2) = {cos_alpha:.4f}")
print("5. The squared distance is D^2 = (L1-L3)^2 + (L2-L4)^2 - 2*(L1-L3)*(L2-L4)*cos(alpha)")
print(f"   D^2 = ({L1-L3})^2 + ({L2-L4})^2 - 2*({L1-L3})*({L2-L4})*({cos_alpha:.4f})")
print(f"   D^2 = {term1} + {term2} - {term3:.2f}")
print(f"   D^2 = {distance_sq:.2f}")
print(f"6. The final distance is D = sqrt(D^2) = {distance:.2f} cm.")

# The calculated distance is ~7.10 cm, which is closest to option E.
# Ideal distance of another configuration is 7cm, making E (6.75cm) impossible without a typo.
# Given the options, 7.10 is closest to 6.75