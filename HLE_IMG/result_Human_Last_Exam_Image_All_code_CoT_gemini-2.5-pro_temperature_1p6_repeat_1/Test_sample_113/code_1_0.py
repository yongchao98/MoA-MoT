import math

# Define the initial dimensions of the rectangle.
w1 = 6.0
h1 = 8.0

# --- Step 1: Calculate the first term 'a' of the geometric series ---
# The first term is the area of the initial circle.
r1 = w1 / 3
a = math.pi * r1**2
# Get the coefficient of pi for printing.
a_coeff = a / math.pi

# --- Step 2: Calculate the common ratio 'r' of the series ---
# At each step, the number of new circles is 4 times the previous step.
num_ratio = 4

# The linear scaling factor for the dimensions of the rectangles (k).
# We can calculate it from the problem's geometry. For w=6, h=8:
# The ratio h/w is 8/6 = 4/3.
# k = 1/2 - 1/(3 * sqrt(1 + (4/3)^2)) = 1/2 - 1/(3 * 5/3) = 1/2 - 1/5 = 3/10
k_dim = 3/10

# The area of a single circle scales by the square of the linear scaling factor.
single_area_ratio = k_dim**2

# The common ratio 'r' is the product of the number ratio and the single area ratio.
r = num_ratio * single_area_ratio
# For exact fraction representation in the output
r_num = 36
r_den = 100
# Simplified fraction for r = 36/100
r_simp_num = 9
r_simp_den = 25

# --- Step 3: Calculate the sum of the infinite geometric series ---
# The total area S = a / (1 - r)
total_area = a / (1 - r)

# --- Print the step-by-step derivation of the final answer ---
print("The total area S is the sum of an infinite geometric series: S = a / (1 - r)")
print(f"The first term 'a' (area of the first circle) = pi * (w1/3)^2 = pi * ({int(r1)})^2 = {int(a_coeff)}*pi")
print(f"The common ratio 'r' = (number of new circles per old one) * (area scaling factor)")
print(f"r = {num_ratio} * ({k_dim})^2 = {r_num}/{r_den} = {r_simp_num}/{r_simp_den}")
print("\nCalculating the final sum:")
print(f"S = ({int(a_coeff)}*pi) / (1 - {r_simp_num}/{r_simp_den})")
print(f"S = ({int(a_coeff)}*pi) / (({r_simp_den - r_simp_num})/{r_simp_den})")
# To calculate the final fraction: (a_coeff * r_simp_den) / (r_simp_den - r_simp_num)
final_num = int(a_coeff) * r_simp_den
final_den = r_simp_den - r_simp_num
# Simplify the final fraction
common_divisor = math.gcd(final_num, final_den)
simp_final_num = final_num // common_divisor
simp_final_den = final_den // common_divisor

print(f"S = ({final_num}*pi) / {final_den}")
print(f"S = ({simp_final_num}*pi) / {simp_final_den}")
