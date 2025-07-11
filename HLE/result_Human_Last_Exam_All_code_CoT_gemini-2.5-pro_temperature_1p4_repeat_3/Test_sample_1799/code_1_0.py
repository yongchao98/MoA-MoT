# The speed of light in a vacuum (meters per second)
c = 299792458

# --- SCENARIO ---
# Imagine a spaceship (Frame S') moving away from a stationary observer on Earth (Frame S).
# The spaceship is moving at a very high velocity, let's say 50% of the speed of light.
v = 0.5 * c

# The spaceship turns on its headlights, emitting a beam of light in its direction of travel.
# An observer on the spaceship will measure the speed of this light beam to be 'c', as expected.

# --- The Question: What speed does the observer on Earth measure? ---

# --- View 1: Classical (Galilean) Relativity ---
# If we only use Postulate 1 and assume classical velocity addition holds true,
# the observer on Earth should see the light moving at its own speed PLUS the speed of the spaceship.
classical_measured_speed = c + v

print("--- Classical (Galilean) Expectation ---")
print("This view assumes only the Principle of Relativity and our everyday intuition for adding speeds.")
print(f"Speed of spaceship (v): {v:,.0f} m/s")
print(f"Speed of light in spaceship's frame (c): {c:,.0f} m/s")
print("Expected speed measured on Earth (u') = c + v")
print(f"u' = {c:,.0f} + {v:,.0f} = {classical_measured_speed:,.0f} m/s\n")


# --- View 2: Einstein's Special Relativity (Using Postulate 2) ---
# Postulate 2 states that the speed of light in a vacuum is the SAME for all inertial observers,
# regardless of the motion of the source.
relativistic_measured_speed = c

print("--- Special Relativity (Einstein's Postulates 1 & 2) ---")
print("This view accepts that the speed of light is constant for all observers.")
print(f"Speed of spaceship (v): {v:,.0f} m/s")
print(f"Speed of light in spaceship's frame (c): {c:,.0f} m/s")
print("Speed measured on Earth according to Postulate 2 is simply c.")
print(f"u' = {relativistic_measured_speed:,.0f} m/s\n")

print("--- Conclusion ---")
print("The two views give different predictions. Experiments have overwhelmingly confirmed that")
print("the prediction from Special Relativity is correct. This shows that the 2nd Postulate")
print("is not superfluous but a necessary principle to correctly describe the universe.")
