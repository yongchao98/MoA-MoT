import math

# --- Problem Parameters ---
R = 1000.0  # Radius of the disk
z0 = (0, 300)  # Starting point
y1 = (0, 0)    # First target point
y2 = (2, 0)    # Second target point

# --- Calculate distances ---
d1 = math.sqrt((z0[0] - y1[0])**2 + (z0[1] - y1[1])**2)
d2 = math.sqrt((z0[0] - y2[0])**2 + (z0[1] - y2[1])**2)

# --- Calculate terms in the formula ---
# Numerator parts
log_R_d1 = math.log(R / d1)
log_R_d2 = math.log(R / d2)
numerator = 4 * (log_R_d1 + log_R_d2)

# Denominator parts
# The term log(15) comes from averaging potentials on the grid.
# The distances from neighbors of y1=(0,0) to y2=(2,0) are 1, 3, sqrt(5), sqrt(5).
# The geometric mean of these distances squared leads to a term log(1*3*5) = log(15).
log_15 = math.log(1 * 3 * 5)
log_R = math.log(R)
denominator = 8 * log_R - log_15

# --- Calculate final probability ---
probability = numerator / denominator

# --- Output the equation and its components ---
print("The formula for the probability P is:")
print(f"P = 4 * (ln(R/d1) + ln(R/d2)) / (8 * ln(R) - ln(15))")
print("\nPlugging in the numbers:")
print(f"R = {R}")
print(f"d1 = ||z0 - y1|| = {d1}")
print(f"d2 = ||z0 - y2|| = {d2:.6f}")
print(f"ln(15) = {log_15:.6f}")
print(f"ln(R) = {log_R:.6f}")
print("\nCalculating the terms:")
print(f"Numerator term 1 (ln(R/d1)): {log_R_d1:.6f}")
print(f"Numerator term 2 (ln(R/d2)): {log_R_d2:.6f}")
print(f"Total Numerator = 4 * ({log_R_d1:.6f} + {log_R_d2:.6f}) = {numerator:.6f}")
print(f"Denominator = 8 * {log_R:.6f} - {log_15:.6f} = {denominator:.6f}")
print("\n--- Final Result ---")
print(f"The probability is P = {numerator:.6f} / {denominator:.6f} = {probability:.6f}")
print(f"\nThe probability rounded to three significant digits is: {probability:.3f}")
