# This script formats and prints the solution to the ordinal arithmetic problem.

# Define string representations for the ordinals for clear output.
w = "ω"
w1 = "ω₁"
w2 = "ω₂"

# The coefficients (alpha values) derived from the step-by-step simplification.
# The final expression is ω₂·ω₁ + ω₁, which corresponds to:
# α₁ = ω₁
# α₂ = 1
# α₃ = 0
# α₄ = 0
a1 = w1
a2 = "1"
a3 = "0"
a4 = "0"

# Print the final expression in the required Cantor Normal Form.
# The request requires printing each number in the final equation.
print("The simplified expression in the form ω₂·α₁ + ω₁·α₂ + ω·α₃ + α₄ is:")
print(f"{w2} · ({a1}) + {w1} · ({a2}) + {w} · ({a3}) + {a4}")

# Also, explicitly state the values of the coefficients for clarity.
print("\nThe values for the ordinals α₁, α₂, α₃, and α₄ are:")
print(f"α₁ = {a1}")
print(f"α₂ = {a2}")
print(f"α₃ = {a3}")
print(f"α₄ = {a4}")