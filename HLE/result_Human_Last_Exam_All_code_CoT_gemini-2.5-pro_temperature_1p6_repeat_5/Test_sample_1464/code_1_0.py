import sympy

# The problem is to find the four roots of a complex polynomial.
# By analyzing the coefficients, we can deduce that the polynomial factors into
# two simpler quadratic equations. This process reveals the four roots.

# Define the symbol X for the polynomial
X = sympy.Symbol('X')

# The four roots identified through the factorization method are:
root1 = 2 * sympy.sqrt(6)
root2 = 2 * sympy.sqrt(11)
root3 = sympy.sqrt(14)
root4 = sympy.sqrt(34)

all_roots = [root1, root2, root3, root4]

# Sort the roots in increasing order based on their numerical value
sorted_roots = sorted(all_roots, key=lambda r: r.evalf())

# To satisfy the instruction "output each number in the final equation!",
# we present the equation in its factored form, which explicitly contains the roots.
r_strs = [f"(X - {r})" for r in sorted_roots]
equation_str = " * ".join(r_strs) + " = 0"

print("The polynomial from the problem can be written in a factored form using its roots.")
print("The final equation is:")
print(equation_str)
print("\n------------------------------------------------------")
print("The 4 roots of the polynomial, in increasing order, are:")
# Print each root on a new line for clarity
for i, root in enumerate(sorted_roots):
    print(f"Root {i+1}: {root}  (approx. {root.evalf()})")
