import sympy

# Define the symbols based on the problem description
d, h, H, rho, t, g = sympy.symbols('d h H rho t g')

# The total mass of the sand
M_s = rho * (sympy.pi * d**2 / 4) * h

# The mass flow rate of the sand
m_dot = M_s / t

# The fall distance for the last grains of sand
L = H - h

# The impact velocity of the last grains of sand
v_impact = sympy.sqrt(2 * g * L)

# The change in weight is estimated by the impact force during the final transient phase
Delta_W = m_dot * v_impact

# The problem requires printing the final equation with each number/symbol separated
# We can use sympy's srepr to get a representation and then format it.
# However, for this specific format, it's easier to construct the string manually
# to match the desired output format of the multiple-choice option.

final_expression_str = "pi * d**2 * h * rho / (4 * t) * sqrt(2 * g * (H - h))"

# Let's format it to look like the option in the text for clarity in the output.
# Option A is: (pi d^2 h rho / 4t) * sqrt(2g(H-h))
# We will print the components of this formula.
# The term for mass flow rate
mass_flow_rate_term = "(pi * d**2 * h * rho) / (4 * t)"
# The term for velocity
velocity_term = "sqrt(2 * g * (H - h))"

# The problem asks to output each number in the final equation. We will format the
# final symbolic expression for printing.
# Let's break down the final expression into parts as requested.
term1 = "pi"
term2 = "d**2"
term3 = "h"
term4 = "rho"
term5 = "4"
term6 = "t"
term7 = "2"
term8 = "g"
term9 = "H"
term10 = "h" # h appears twice

print(f"The change in weight, Delta W, can be estimated by the final impact force spike.")
print(f"The expression for the change in weight is:")
# Manually constructing the print string to match the format of option A
# It's a bit ambiguous what "output each number" means for a symbolic expression,
# so we will print the formula in a readable way that shows all its components.
print(f"Delta_W = ( {term1} * {term2} * {term3} * {term4} ) / ( {term5} * {term6} ) * sqrt( {term7} * {term8} * ( {term9} - {term10} ) )")

# For direct comparison with the options, let's print the simplified formula string.
# Option A: A. \frac{\pi d^2 h \rho}{4t}\sqrt{2g(H-h)}
# We will print the python equivalent of this.
print("\nOr as a single line expression:")
print("Delta_W = (pi * d**2 * h * rho / (4 * t)) * sqrt(2 * g * (H - h))")
