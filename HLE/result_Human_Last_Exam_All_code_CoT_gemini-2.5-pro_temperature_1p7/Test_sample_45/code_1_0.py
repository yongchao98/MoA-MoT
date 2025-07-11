# Step 1: Define the relationship between focal length (f) and angular speed (omega).
# The surface of a rotating liquid is a paraboloid z = omega^2 * x^2 / (2 * g).
# A parabolic mirror has the equation x^2 = 4 * f * z, so z = x^2 / (4 * f).
# Comparing these gives 1 / (4 * f) = omega^2 / (2 * g), so f = g / (2 * omega^2).
# This means f is proportional to 1 / omega^2, or f ~ (omega^2)^(-1).

# Step 2: Define the relationship between angular speed (omega) and time (t).
# The system is driven by a constant power source (P).
# Power is the rate of change of rotational kinetic energy (Ek): P = d(Ek)/dt.
# Ek = 0.5 * I * omega^2, where I is the moment of inertia.
# P = d(0.5 * I * omega^2)/dt = I * omega * d(omega)/dt.
# Separating variables: P * dt = I * omega * d(omega).
# Integrating from t=0 (omega=0) gives P * t = 0.5 * I * omega^2.
# So, omega^2 = (2 * P / I) * t.
# This means omega^2 is proportional to t.

# Step 3: Combine the relationships.
# We have f ~ (omega^2)^(-1) and omega^2 ~ t.
# Substituting the second into the first gives f ~ t^(-1).

# Step 4: Determine the value of n.
# The problem states f is proportional to t^n.
# From our derivation, the exponent n is -1.
n = -1

# Print the final answer.
print("The focal length f is proportional to t^n.")
print(f"The value of n is: {n}")