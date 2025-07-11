# Define the integer frequencies for the Lissajous curve (cos(nx*t), sin(ny*t))
nx = 9
ny = 5

# The number of self-intersection points for a Lissajous curve with coprime
# frequencies p and q is given by the formula (p - 1) * (q - 1) / 2.
# In this case, nx=9 and ny=5 are coprime, so we use p=9 and q=5.

# Calculate the terms of the formula
p_minus_1 = nx - 1
q_minus_1 = ny - 1
numerator = p_minus_1 * q_minus_1
denominator = 2
result = numerator // denominator

# Print the step-by-step evaluation of the formula
print(f"The number of self-intersection points for the curve (cos({nx}t), sin({ny}t)) is calculated as follows:")
print(f"Formula: (p - 1) * (q - 1) / 2")
print(f"Substituting p={nx} and q={ny}:")
print(f"({nx} - 1) * ({ny} - 1) / {denominator} = {p_minus_1} * {q_minus_1} / {denominator} = {numerator} / {denominator} = {result}")
