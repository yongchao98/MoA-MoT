# Define the string representations for the ordinals
w2 = "ω₂"
w1 = "ω₁"
w = "ω"

# Define the coefficients based on the derivation
# The final expression is w2*w1 + w1*1
# Which maps to w2*a1 + w1*a2 + w*a3 + a4
a1 = "ω₁"
a2 = "1"
a3 = "0"
a4 = "0"

# Print the final expression in the requested format
# The full form is ω₂⋅α₁ + ω₁⋅α₂ + ω⋅α₃ + α₄
print(f"The simplified expression is:")
print(f"{w2} * {a1} + {w1} * {a2} + {w} * {a3} + {a4}")