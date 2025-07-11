import sympy as sp

# Define symbols for our physical quantities.
# t: time
# g: acceleration due to gravity
# P: power from the source
# I: moment of inertia of the liquid
t, g, P, I = sp.symbols('t g P I', positive=True, real=True)
omega_sq = sp.Symbol('omega_sq', positive=True, real=True) # Represents ω²

# --- Derivation Steps ---
print("--- Derivation Steps ---")

# Step 1: Find focal length f in terms of angular velocity ω.
# The surface of a rotating liquid is a parabola z = ω²r²/(2g).
# The focal length 'f' for a standard parabola y = x²/(4f) leads to the relation f = g/(2ω²).
f_from_omega_sq = g / (2 * omega_sq)
print("1. The focal length f is related to the angular velocity squared ω² by the equation:")
print(f"f = g / (2 * ω²)")
print("The number in this equation is 2.")
print("-" * 30)

# Step 2: Find angular velocity ω in terms of time t.
# A constant power source P means the kinetic energy K increases linearly with time: K = P*t (starting from rest).
# The rotational kinetic energy is K = (1/2)*I*ω².
# Equating these gives (1/2)*I*ω² = P*t, which we solve for ω².
omega_sq_from_t = (2 * P * t) / I
print("2. With a constant power source, ω² is related to time t by the equation:")
print(f"ω² = (2 * P * t) / I")
print("The number in this equation is 2.")
print("-" * 30)

# Step 3: Combine the relations to find f as a function of t.
# We substitute the expression for ω² from step 2 into the equation from step 1.
f_of_t = f_from_omega_sq.subs(omega_sq, omega_sq_from_t)
f_of_t_simplified = sp.simplify(f_of_t) # Result is g*I/(4*P*t)
print("3. Combining these two relationships gives the focal length f as a function of time t:")
print(f"f(t) = {f_of_t_simplified}")
print("\nThis can be rewritten in the form f(t) = C * tⁿ as:")
print(f"f(t) = (g*I / (4*P)) * t**(-1)")
print("The numbers in this final expression are 4 in the denominator and the exponent -1.")
print("-" * 30)


# Step 4: Determine the exponent n.
# The problem asks for the value of n, where f is proportional to tⁿ (f ∝ tⁿ).
# From our result, f(t) is proportional to t⁻¹, so we can directly identify n.
n = -1

print("\n--- Final Answer ---")
print("The problem asks to compute n for the relation f ∝ tⁿ.")
print("Based on the derivation, the final equation for the exponent is:")
print(f"n = {n}")