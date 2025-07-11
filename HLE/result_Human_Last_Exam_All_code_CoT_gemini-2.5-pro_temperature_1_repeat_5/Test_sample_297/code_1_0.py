import sympy
from sympy.logic.boolalg import And, Or, to_dnf

# This script calculates the number of multiplications in the fully expanded
# expression for the sum bit s2 of a 3-bit binary adder.

# 1. Define symbolic boolean variables for the inputs.
# The complement of a variable 'a' will be represented by '~a'.
a2, a1, a0, b2, b1, b0 = sympy.symbols('a2 a1 a0 b2 b1 b0', cls=sympy.logic.boolalg.BooleanSymbol)

# 2. Define the expression for the first carry bit, c1.
c1 = And(a0, b0)

# 3. Define the expression for the second carry bit, c2.
c2 = Or(And(a1, b1), And(a1, c1), And(b1, c1))

# 4. Define the expression for the sum bit, s2, using the XOR operator (^).
s2 = a2 ^ b2 ^ c2

# 5. Substitute c1 into c2, and then the resulting c2 into s2 to get the
#    full expression for s2 in terms of only the primary input bits.
s2_full_expression = s2.subs(c2, c2.subs(c1, And(a0, b0)))

# 6. Convert the expression to its Disjunctive Normal Form (DNF), which is
#    the "fully expanded" sum-of-products form.
s2_dnf = to_dnf(s2_full_expression, simplify=True)

# 7. Count the multiplication operations. In a sum-of-products expression,
#    OR is "addition" and AND is "multiplication". A product term with 'n'
#    variables has 'n-1' multiplications.
total_multiplications = 0
equation_terms = []

def format_term(term):
    """Formats a sympy boolean expression for readable output (e.g., a2b2'a1)."""
    # Replace sympy's Not symbol (~) with a prime (') and remove the AND operator (&).
    s = str(term).replace('~', "'").replace(' & ', '')
    return s

# The result of to_dnf can be an Or, an And, or a single literal.
# We create a list of terms to handle all cases uniformly.
if isinstance(s2_dnf, sympy.logic.boolalg.Or):
    terms = s2_dnf.args
else:
    terms = [s2_dnf]

for term in terms:
    # A term in a DNF can be an 'And' object or a single variable.
    if isinstance(term, sympy.logic.boolalg.And):
        num_variables = len(term.args)
        # A product of n variables has n-1 multiplications.
        if num_variables > 1:
            total_multiplications += (num_variables - 1)
    # A single literal term has 0 multiplications.

    # Format the term to display in the final equation.
    equation_terms.append(format_term(term))

# 8. Print the final results, including the equation and the total count.
# The equation shows each digit (0, 1, 2) in the variable names.
final_equation_str = "s2 = " + " + ".join(sorted(equation_terms))
print("The fully expanded expression for s2 is:")
print(final_equation_str)

print("\nTotal number of multiplication operations in the expression:")
print(total_multiplications)
<<<52>>>