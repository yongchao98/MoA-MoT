# The following code prints the final result of the ordinal simplification.
# Unicode characters are used for mathematical symbols for clarity.
# ω -> \u03C9
# ₁ -> \u2081
# ₂ -> \u2082
# ⋅ -> \u22C5
# α -> \u03B1

# The simplified expression is ω₂ ⋅ ω₁ + ω₁ ⋅ 1.
# We present it in the full requested form ω₂⋅α₁ + ω₁⋅α₂ + ω⋅α₃ + α₄
# by making the coefficients explicit.

final_expression = "\u03C9\u2082 \u22C5 (\u03C9\u2081) + \u03C9\u2081 \u22C5 (1) + \u03C9 \u22C5 (0) + 0"

print("The original expression simplifies to the following form:")
print(final_expression)
print("Where the coefficients are:")
print("  \u03B1\u2081 = \u03C9\u2081")
print("  \u03B1\u2082 = 1")
print("  \u03B1\u2083 = 0")
print("  \u03B1\u2084 = 0")
