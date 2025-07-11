import math

# Use fractions to emulate the decimal-based frac type arithmetic
from fractions import Fraction

# 1. Define constants from the problem
luminosity = 2e9
mass_light_ratio = 3
velocity = 200
radius = 10

# 2. Calculate Luminous Mass
m_lum = luminosity * mass_light_ratio

# 3. Calculate Total Mass
# Numerator v^2 * r
v2r = velocity**2 * radius
# Denominator G, approximated as 13/3 * 10^-6
g_approx = Fraction(13, 3) * 1e-6
# Total Mass
m_tot = v2r / g_approx

# 4. Calculate Dark Matter Percentage
dark_matter_ratio = 1 - (m_lum / m_tot)
dark_matter_percentage = dark_matter_ratio * 100

# 5. Format the numbers for the equation string
# Luminous Mass: 6e9 -> 6.00e+09
# Total Mass: M_tot value -> formatted string
# Percentage: Rounded to one decimal place.
equation_str = (
    f"Dark Matter Percentage = (1 - {m_lum:.2e} / {m_tot:.2e}) * 100"
    f" = {dark_matter_percentage:.1f}%"
)

# Print each number in the final equation, as requested.
print(equation_str)

# 6. Calculate the memory usage (z)
# 4 frac variables * 6D/frac + 1 int variable * 5D/int
z = 4 * 6 + 1 * 5
p = f"{dark_matter_percentage:.1f}"

# The final answer must be in the format p:z
# print(f"Final Answer in p:z format is {p}:{z}")