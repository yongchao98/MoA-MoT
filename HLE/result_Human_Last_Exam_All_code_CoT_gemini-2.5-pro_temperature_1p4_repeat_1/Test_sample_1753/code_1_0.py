import math

# We need to solve the equation for the total arc length:
# 2 * Integral from t_a to pi/2 of (3/2)sin(2t) dt = 3/2
# This simplifies to: Integral from t_a to pi/2 of sin(2t) dt = 1/2

# 1. Set up the equation from the integral
# Integral(sin(2t)) = -1/2 * cos(2t)
# Evaluating from t_a to pi/2:
# (-1/2 * cos(pi)) - (-1/2 * cos(2*t_a)) = 1/2
# (-1/2 * -1) + (1/2 * cos(2*t_a)) = 1/2
val_pi_term = 1/2
eq_rhs = 1/2
print(f"The equation derived from the integral is: {val_pi_term} + (1/2)*cos(2*t_a) = {eq_rhs}")

# 2. Solve for cos(2*t_a)
cos_2ta = (eq_rhs - val_pi_term) * 2
print(f"This simplifies to cos(2*t_a) = {cos_2ta}")

# 3. Use the double angle identity cos(2t) = 2*cos(t)^2 - 1
# 2*cos(t_a)^2 - 1 = 0
eq_rhs_cos_sq = 1
print(f"Using the identity cos(2t) = 2*cos(t)^2 - 1, we get: 2*cos(t_a)^2 = {eq_rhs_cos_sq}")

# 4. Solve for cos(t_a)
cos_ta_sq = eq_rhs_cos_sq / 2
print(f"Therefore, cos(t_a)^2 = {cos_ta_sq}")
# Since t_a is in the first quadrant, cos(t_a) is positive.
cos_ta = math.sqrt(cos_ta_sq)
print(f"Since t_a is in the first quadrant, cos(t_a) = 1/sqrt(2) ≈ {cos_ta:.7f}")

# 5. Solve for 'a' using the relation a = cos(t_a)^3
a = cos_ta ** 3
print(f"From the relation a = cos(t_a)^3, we find the value of a:")
print(f"a = (1/sqrt(2))^3 = 1/(2*sqrt(2))")
print(f"The numerical value of a is approximately {a:.7f}")
final_a = a
<<<0.3535533905932738>>>