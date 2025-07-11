# The purpose of this script is to generate and print the two additional inequalities
# required for the MILP formulation.

# The parameters and variables are represented as string symbols.
y = "y"
x = "x"
a = "a"
b = "b"
l = "l"
u = "u"

# The derived inequalities are:
# 1. y >= l*a - l*b + l
# 2. y >= x - u*a - u*b
# These are already in a distributed form. We will construct them using f-strings.

inequality_1 = f"{y} >= {l}*{a} - {l}*{b} + {l}"
inequality_2 = f"{y} >= {x} - {u}*{a} - {u}*{b}"

# The final output is formatted as a single line with the two inequalities separated by a comma.
print(f"{inequality_1}, {inequality_2}")
