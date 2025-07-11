import math

print("The problem statement is ambiguous, leading to two possible values for 'a'.")
print("Here are the derivations for both interpretations.\n")

# --- Interpretation 1: The arc is a single continuous path from x=0 to x=a. ---
print("--- Interpretation 1: Arc is a continuous path from x=0 to x=a ---")

# The relationship between the arc length L and the endpoint 'a' is derived as:
# L = (3/2) * cos^2(t_a), where 'a' is related to the parameter t_a by a = cos^3(t_a).
L_given_1 = 3/2
coefficient_1 = 3/2

# We are given L = 3/2. We set the equations equal to solve for 'a'.
# Equation: L_given = coefficient * cos^2(t_a)
# -> 3/2 = (3/2) * cos^2(t_a)
cos_t_a_squared = L_given_1 / coefficient_1
# cos(t_a) must be non-negative for x to be non-negative.
cos_t_a = math.sqrt(cos_t_a_squared)
# The value for a is then calculated from a = cos^3(t_a).
a1 = cos_t_a**3

print(f"Given arc length L = {L_given_1}")
print(f"The derived length formula is: L = {coefficient_1} * cos^2(t_a)")
print(f"Equating them gives: {L_given_1} = {coefficient_1} * cos^2(t_a)")
print(f"This simplifies to: cos^2(t_a) = {L_given_1 / coefficient_1}")
print(f"So, cos(t_a) = sqrt({cos_t_a_squared}) = {cos_t_a}")
print(f"From a = cos^3(t_a), we get a = {cos_t_a}^3 = {a1}")
print(f"Value for 'a' under Interpretation 1: {a1}\n")


# --- Interpretation 2: The arc is the set of all points where 0 <= x <= a. ---
print("--- Interpretation 2: Arc is the locus of points where 0 <= x <= a ---")

# The total length L of the curve portion where 0 <= x <= a is derived as:
# L = 3 * a^(2/3)
L_given_2 = 3/2
coefficient_2 = 3

# We are given L = 3/2. We set the equations equal to solve for 'a'.
# Equation: L_given = coefficient * a^(2/3)
# -> 3/2 = 3 * a^(2/3)
a_pow_2_3 = L_given_2 / coefficient_2
# The value for a is calculated by solving a^(2/3) = 1/2.
a2 = a_pow_2_3**(3/2)

print(f"Given arc length L = {L_given_2}")
print(f"The derived length formula is: L = {coefficient_2} * a^(2/3)")
print(f"Equating them gives: {L_given_2} = {coefficient_2} * a^(2/3)")
print(f"This simplifies to: a^(2/3) = {L_given_2} / {coefficient_2} = {a_pow_2_3}")
print(f"Solving for a: a = ({a_pow_2_3})^(3/2)")
print(f"Value for 'a' under Interpretation 2: {a2} (which is sqrt(2)/4)")
