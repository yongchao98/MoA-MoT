import math

# Time t
t = math.pi / 8

# Value of z1(t) at t=pi/8
# z1(t) = cos(t)
z1_at_t = math.cos(t)

# Value of y(t) at t=pi/8
# y(t) = tan(2t) * sqrt(sec(2t))
tan_2t = math.tan(2 * t)
sec_2t = 1 / math.cos(2 * t)
y_at_t = tan_2t * math.sqrt(sec_2t)

# Calculate the final expression (z1(t) / y(t))^2
result = (z1_at_t / y_at_t)**2

# Print the components and the final equation with numerical values
print("The problem requires calculating (z1(π/8) / y(π/8))^2.")
print(f"The value of the numerator term is z1(π/8) = cos(π/8) ≈ {z1_at_t:.6f}")
print(f"The value of the denominator term is y(π/8) = tan(π/4) * sqrt(sec(π/4)) ≈ {y_at_t:.6f}")
print("\nThe final equation with these values is:")
print(f"({z1_at_t:.6f} / {y_at_t:.6f})^2 = {result:.6f}")
print("\nThe exact analytical result is (√2 + 1)/4.")
