from fractions import Fraction

# Step 1: Define initial values from the problem analysis
# v_a5 is the valuation of the leading coefficient of f(x), which is 2.
v_a5 = Fraction(1, 1)

# Step 2: Calculate the valuation of the derivative at a root alpha.
# From the analysis: v(f'(alpha)) = v(2) + v(alpha^2) + v(5*alpha^2 + 3)
v_alpha = Fraction(-1, 5)
v_alpha_sq = 2 * v_alpha
v_5_alpha_sq_plus_3 = v_alpha_sq # As v(5*alpha_sq) = -2/5 < v(3)=0
v_f_prime = v_a5 + v_alpha_sq + v_5_alpha_sq_plus_3

# Step 3: Calculate V, the sum of valuation differences.
# The formula is v(f'(alpha)) = v(a_n) + V
V = v_f_prime - v_a5

# Step 4: Calculate the thickness.
# The formula is thickness = -1/2 * V
factor = Fraction(-1, 2)
thickness = factor * V

# Output the steps of the final calculation
print("The thickness of the double point is calculated as follows:")
print("thickness = (-1/2) * (v(f'(alpha)) - v(a_5))")
print(f"where v(a_5) = v(2) = {v_a5}")
print(f"and v(f'(alpha)) = v(2) + v(alpha^2) + v(5*alpha^2+3) = {v_a5} + {v_alpha_sq} + {v_5_alpha_sq_plus_3} = {v_f_prime}")
print(f"So, the calculation is:")
print(f"thickness = {factor} * ({v_f_prime} - {v_a5}) = {thickness}")
