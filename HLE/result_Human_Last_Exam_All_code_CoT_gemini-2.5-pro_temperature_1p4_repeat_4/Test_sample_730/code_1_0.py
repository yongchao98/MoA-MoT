# The components of the calculation are determined by the problem's data and formulas.
# Density (ρ) = 0.9 = 9/10
density_num = 9
density_den = 10

# Volume constant = 4/3
volume_const_num = 4
volume_const_den = 3

# π is approximated as 3, which has an error of about 4.5% (< 10%).
pi_approx = 3

# Radius (r) = 0.5 = 1/2. The formula uses r^3, which is (1/2)^3 = 1/8.
radius_cubed_num = 1
radius_cubed_den = 8

# The final instruction for the parrot is the full calculation.
# We print each number involved in the equation.
print("Here is the calculation to estimate the mass:")
print(f"Mass = {density_num}/{density_den} * {volume_const_num}/{volume_const_den} * {pi_approx} * {radius_cubed_num}/{radius_cubed_den}")