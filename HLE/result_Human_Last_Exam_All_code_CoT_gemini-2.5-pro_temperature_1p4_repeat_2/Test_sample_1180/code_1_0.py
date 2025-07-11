import sympy

# We are given the valuation v(2) = 1.
v_2 = 1

# From analyzing the singularity of the curve's model, we derived a system of
# linear equations for the valuations of the local coordinates X and Z1.
# Let v(X) be x and v(Z1) be y.
# The equations are:
# 1) 2*y - x = v(2)
# 2) y - 2*x = 0
x, y = sympy.symbols('x y')
solution = sympy.solve([
    2*y - x - v_2,
    y - 2*x
], (x, y))

v_X = solution[x]
v_Z1 = solution[y]

# The thickness of the double point is the common valuation of the
# principal terms defining the singularity.
# Thickness = v(2*X) = v(2) + v(X)
thickness = v_2 + v_X

# We output the equation used to calculate the thickness
# along with its final value.
print(f"The thickness is calculated as v(2) + v(X).")
print(f"With v(2) = {v_2} and the calculated v(X) = {v_X}, the equation is:")
# sympy.Rational automatically finds the common denominator for printing.
print(f"{v_2} + {v_X} = {thickness}")
