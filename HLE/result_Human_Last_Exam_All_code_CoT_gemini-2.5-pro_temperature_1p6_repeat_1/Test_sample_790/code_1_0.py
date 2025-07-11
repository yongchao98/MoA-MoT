import math

# Titan Computer Simulation for Coconut Problem

# --- Step 1: Define Constants and Approximations ---
# All numbers must be represented as fractions a/b where a, b are in [0, 31].

# Rock radius (r) = 0.5 cm -> 1/2
r_num, r_den = 1, 2
# Rock density (rho) = 0.9 kg/cm^3 -> 9/10
rho_num, rho_den = 9, 10

# Gravitational acceleration (g) approx 9.8 m/s^2.
# We approximate g as 10/1 to simplify calculations and avoid overflows.
# 9.8 = 49/5 (49>31). 29/3=9.67 is also possible but creates overflows later.
# 10/1 has a ~2% error but makes the calculation feasible.
g_num, g_den = 10, 1

# Mathematical constant pi approx 3.14159.
# We approximate pi as 3/1. Using 22/7 leads to intermediate overflows.
# This is a significant approximation (~4.5% error) but necessary.
pi_num, pi_den = 3, 1

# Launch angle is 45 degrees. We need cos(45) = sqrt(2)/2.
# From the problem, we can use a high-precision approximation for sqrt(2) ~ 24/17.
# This gives cos(45) ~ (24/17) / 2 = 12/17.
cos45_num, cos45_den = 12, 17

print("Titan Computer Approximations:")
print(f"g = {g_num}/{g_den}")
print(f"pi = {pi_num}/{pi_den}")
print(f"cos(45) = {cos45_num}/{cos45_den}\n")


# --- Step 2: Calculate Mass (m) of the Rock ---
# m = Volume * density = (4/3) * pi * r^3 * rho
# m = (4/3) * (3/1) * (1/2)^3 * (9/10)
# Let's do this step-by-step to respect Titan's constraints.
# term1 = (4/3) * (3/1) = 4/1
# term2 = (1/2)^3 = 1/8
# term3 = (9/10)
# m = (4/1) * (1/8) * (9/10)
# m = (1/2) * (9/10) = 9/20
m_num, m_den = 9, 20
print(f"Calculated mass (m) = {m_num}/{m_den} kg\n")


# --- Step 3: Calculate the Force (F) ---
# The physics formula for the required force is F = 2 * m * g / cos(45)
print("Calculating force F = (2/1) * m * g * (1/cos(45))...")

# Combine terms sequentially, simplifying at each step.
# F = (2/1) * (9/20) * (10/1) * (17/12)
# F_step1 = (2/1) * (9/20) = 18/20 = 9/10
# F_step2 = (9/10) * (10/1) = 90/10 = 9/1
# F = (9/1) * (17/12)

F_interim_num, F_interim_den = 9, 1
op_num, op_den = 17, 12
print(f"Intermediate calculation: F = ({F_interim_num}/{F_interim_den}) * ({op_num}/{op_den})")

# Now we must perform the multiplication (9/1) * (17/12).
# The new numerator would be 9 * 17 = 153.
# This exceeds the 5-bit limit of 31. The operation is invalid on Titan.
print(f"Invalid operation: Numerator {F_interim_num}*{op_num} = {F_interim_num*op_num} which is > 31.")

# Per the Titan rules, we can substitute a fraction with a close approximation
# to make the calculation possible.
# Value of 17/12 is approx 1.4166.
# A good approximation with small integers is 13/9, which is approx 1.4444.
# This is a valid strategy to maintain computability.
approx_num, approx_den = 13, 9
print(f"Approximating {op_num}/{op_den} with {approx_num}/{approx_den} to proceed.\n")

# --- Step 4: Final Calculation ---
# F = (9/1) * (13/9)
# The 9s cancel out, resulting in 13/1.
final_F_num = 13
final_F_den = 1

print("Final Result:")
# As requested, outputting each number in the final equation.
print(f"({F_interim_num} / {F_interim_den}) * ({approx_num} / {approx_den}) = {final_F_num} / {final_F_den}")
print(f"The required force is {final_F_num}/{final_F_den} N.")

# --- Bonus: Error Calculation ---
# "True" value using floating point math for comparison
true_pi = math.pi
true_g = 9.8
true_m = (3/20) * true_pi # approx 0.4712
true_cos45 = math.cos(math.radians(45)) # approx 0.7071
true_F = (2 * true_m * true_g) / true_cos45 # approx 13.062

our_F = final_F_num / final_F_den
error = abs(our_F - true_F)
print(f"\nThe calculated force is {our_F:.3f} N.")
print(f"A more precise calculation gives {true_F:.3f} N.")
print(f"The absolute error is {error:.3f}.")
