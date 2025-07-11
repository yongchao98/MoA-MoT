import sympy

# Define the symbols for the side lengths of the triangle ABC.
# a represents the length of side BC (opposite to vertex A).
# b represents the length of side CA (opposite to vertex B).
# c represents the length of side AB (opposite to vertex C).
a, b, c = sympy.symbols('a b c')

# The derivation using geometric properties (Ptolemy's Theorem and the Incenter-Excenter Lemma)
# provides the expression for the ratio BM/MI in terms of the side lengths a, b, and c.

# The numerator of the expression is the sum of the sides adjacent to angle B.
numerator = a + c

# The denominator is the side opposite to angle B.
denominator = b

# Now, we print the final equation for the ratio.
# The instruction to "output each number in the final equation" is interpreted as
# displaying the complete symbolic formula clearly, showing all its components.
print(f"The ratio BM / MI is given by the following equation:")

# The f-string formats the output by converting the sympy symbolic expressions to strings.
# This results in a clear representation of the final equation.
print(f"BM / MI = ({numerator}) / ({denominator})")