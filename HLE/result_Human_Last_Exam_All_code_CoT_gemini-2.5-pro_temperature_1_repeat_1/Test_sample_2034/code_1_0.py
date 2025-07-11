# The determined values for the coefficients c1 and c2.
c1 = 1
c2 = 0

# Print the determined values for the user.
print(f"The determined value for c1 is: {c1}")
print(f"The determined value for c2 is: {c2}")

# Print the final form of the inequality with the specific numbers.
# We use LaTeX formatting for clarity as in the problem description.
print("\nThe final inequality is:")
final_equation = f"-(K \\vec{{\\alpha}}^{{\\mathcal{{D}}-i}})_i \\leq (1 + {c1}\\beta)\\alpha^{{\\mathcal{{D}}}}_i - (1 + {c2}\\beta)(K \\vec{{\\alpha}}^{{\\mathcal{{D}}}}) _i + o(\\beta)"
print(final_equation)

# Print the simplified version of the inequality.
print("\nWhich simplifies to:")
simplified_equation = f"-(K \\vec{{\\alpha}}^{{\\mathcal{{D}}-i}})_i \\leq (1 + \\beta)\\alpha^{{\\mathcal{{D}}}}_i - (K \\vec{{\\alpha}}^{{\\mathcal{{D}}}}) _i + o(\\beta)"
print(simplified_equation)