import sympy
from sympy import pi, sqrt

# Define the symbols
d, h, H, rho, t, g = sympy.symbols('d h H rho t g')

# Mass of the sand
M_sand = rho * (pi * d**2 / 4) * h

# Mass flow rate
m_dot = M_sand / t
m_dot_expr = (pi * d**2 * h * rho) / (4 * t)

# Maximum impact velocity (falling from height H)
v_impact_max = sqrt(2 * g * H)

# Change in weight is estimated by the largest possible impact force
Delta_W = m_dot * v_impact_max
Delta_W_expr = m_dot_expr * v_impact_max

# Print the final expression
# sympy.pretty_print is not available, so we construct the string manually
# For \frac{\pi d^2 h \rho}{4t}\sqrt{2g(H-h)}, the latex is:
# \frac{\pi d^2 h \rho}{4t}\sqrt{2gH}
# Let's print the components to form the final expression

print("The change in weight, Delta W, is estimated by the largest possible impact force.")
print("This occurs when the sand falls the maximum height H.")
print("\nStep 1: The mass flow rate (m_dot) is the total mass of sand divided by the time it takes to fall.")
print(f"m_dot = (pi * d^2 * h * rho) / (4 * t)")
print("\nStep 2: The maximum impact velocity (v_max) occurs after falling a distance H.")
print(f"v_max = sqrt(2 * g * H)")
print("\nStep 3: The estimated weight change is the product of the mass flow rate and the maximum velocity.")
print(f"Delta W = m_dot * v_max")
print("\nFinal Expression:")
# We want to print the formula like in the options
# Option D is (pi * d**2 * h * rho)/(4*t) * sqrt(2*g*H)
# We can print it by concatenating the string representation of the sympy objects
print(f"Delta W = ({sympy.pretty(m_dot_expr, use_unicode=False)}) * ({sympy.pretty(v_impact_max, use_unicode=False)})")
print(f"Delta W = (pi*d**2*h*rho)/(4*t) * sqrt(2*g*H)")
print("\nThis means the hourglass is heavier while running. The final equation is:")
print(f"Delta W = (pi * d**2 * h * rho / (4 * t)) * sqrt(2 * g * H)")
# To match the output format of the choices
print(f"Final equation: Delta W = (pi * d^2 * h * rho / (4*t)) * sqrt(2*g*H)")
# Let's print the numbers in the final equation: pi, 2, 4
print("The numbers in the equation are:")
print(f"Numerator constant: {1}")
print(f"Denominator constant: {4}")
print(f"Constant in sqrt: {2}")