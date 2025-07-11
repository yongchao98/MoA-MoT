# This script calculates the exponent n in the relationship f ∝ t^n for a liquid-mirror telescope
# driven by a constant power source. The calculation is broken down into steps derived from physics principles.

# Step 1: Find the exponent 'p' in the relationship f ∝ (ω^2)^p
# The focal length 'f' of a parabolic mirror created by a rotating liquid is given by f = g / (2 * ω^2),
# where 'g' is the acceleration due to gravity and 'ω' is the angular speed.
# From this equation, we can see that the focal length is inversely proportional to the square of the angular speed.
# This relationship can be written as f ∝ 1 / ω^2, which is equivalent to f ∝ (ω^2)^(-1).
# So, the exponent p is -1.
p = -1

# Step 2: Find the exponent 'q' in the relationship ω^2 ∝ t^q
# The system is driven by a constant power source P. Power is the rate of change of kinetic energy K, so P = dK/dt.
# The rotational kinetic energy is K = (1/2) * I * ω^2, where 'I' is the constant moment of inertia.
# The equation for power becomes d/dt( (1/2) * I * ω^2 ) = P.
# Since I and P are constants, we can integrate with respect to time 't'.
# Assuming the rotation starts from rest (ω=0 at t=0), the integrated equation is (1/2) * I * ω^2 = P * t.
# This shows that ω^2 is directly proportional to time t, or ω^2 ∝ t^1.
# So, the exponent q is 1.
q = 1

# Step 3: Combine the exponents to find n in the relation f ∝ t^n
# We have f ∝ (ω^2)^p and ω^2 ∝ t^q.
# By substituting the second relation into the first, we get f ∝ (t^q)^p = t^(p*q).
# Therefore, the final exponent is n = p * q.
n = p * q

print("To find the exponent 'n' in the relation f ∝ t^n, we combine two relationships:")
print(f"1. The focal length vs. angular speed relation is: f ∝ (ω^2)^p")
print(f"   The value of the exponent p is: {p}")
print("")
print(f"2. The angular speed vs. time relation is: ω^2 ∝ t^q")
print(f"   The value of the exponent q is: {q}")
print("")
print("3. The final relationship is f ∝ t^n, where n = p * q.")
print("The final equation for n is based on the numbers derived above:")
print(f"   n = ({p}) * ({q})")
print(f"   n = {n}")
