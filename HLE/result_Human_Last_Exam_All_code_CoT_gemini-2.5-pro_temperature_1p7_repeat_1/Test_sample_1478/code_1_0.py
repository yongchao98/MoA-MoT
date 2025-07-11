# This program simulates the Titan 6-bit architecture to determine
# if the specified physics calculation is feasible by checking for overflows.

# --- Step 1: Define the Problem and Constants ---
# The physics formula for the force F is approximated as:
# F ≈ (G * rho * (4/3) * pi * R_planet^3 * m) / d^2
# We must use fractional representations for constants where numerators
# and denominators are within the 6-bit limit (<= 63).

print("Formula for Force (F) coefficients:")
print("F_coeff = G_c * rho_c * (4/3) * pi_c * R_c^3 * m_c / d_c^2")
print("")
print("Approximations for coefficients:")
print("G_c = 20/3")
print("rho_c = 12/1")
print("4/3 = 4/3")
print("pi_c = 22/7")
print("R_c^3 = 8/1")
print("m_c = 50/1")
print("d_c^2 = 1/1")
print("")

# --- Step 2: Simulate the Titan Calculation Step-by-Step ---
# The Titan computer must perform multiplications sequentially. Any operation
# that results in a numerator or denominator > 63 is an overflow and fails.

print("--- Simulating Titan Operations ---")

# Let's start by calculating the term for the planet's mass properties.
# A reasonable first step is multiplying the density by the 4/3 volume factor.
print("MOV AX, 12/1  # Load density mantissa (rho_c)")
print("MUL AX, 4/3   # Multiply by volume factor")

# Calculation: 12/1 * 4/3 = 48/3.
# We can simplify this fraction by dividing both by their GCD(3), which results in 16/1.
# Both 16 and 1 are <= 63, so this step is valid.
print("# AX holds 48/3, which reduces to 16/1. This intermediate step is valid.")
print("# Current value in AX: 16/1")
print("")

# Now, we must multiply this result by another term. Let's try the cube of the
# radius's coefficient, which is 8.
print("MUL AX, 8/1   # Attempt to multiply by the radius-cubed mantissa (R_c^3)")

# Calculation: 16/1 * 8/1
numerator = 16 * 8
denominator = 1 * 1
print(f"# This operation would result in: {numerator}/{denominator}")

# Here, we check for the overflow condition.
MAX_VAL = 63
if numerator > MAX_VAL or denominator > MAX_VAL:
    print(f"# FAILURE: The resulting numerator, {numerator}, is greater than {MAX_VAL}.")
    print("# This operation violates the 6-bit constraint.")
else:
    print("# This step is valid.")

# --- Step 3: Conclusion ---
# We have shown that a necessary multiplication step results in an overflow.
# Exploring other multiplication orders yields similar failures:
# - m_c * R_c^3   = 50 * 8   = 400 (>63) -> FAIL
# - rho_c * R_c^3 = 12 * 8   = 96 (>63)  -> FAIL
# - (Result of first step: 16) * pi_c = 16 * 22/7 = 352/7 (>63) -> FAIL

print("")
print("# No matter the order of operations, the product of the given coefficients")
print("# cannot be calculated without intermediate numerators exceeding 63.")
print("# Therefore, the task is impossible to perform on the Titan computer.")

print("<<<N0>>>")