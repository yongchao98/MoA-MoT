import numpy as np

# The value of alpha is determined by the logic above.
# We can express sqrt(2) as a number.
alpha = np.sqrt(2)

# The final answer must be a numerical value.
print(f"The determined value of alpha is the square root of 2.")
print(f"sqrt(2) = {alpha}")
# Final equation: The growth rate of c_n as n to infinity can be written c_n = Theta( (sqrt(2))^n ).
# The base of the exponent is alpha.
# Here we print the numbers required for the equation: c_n, n, alpha
print(f"c_n = Theta(({alpha})^n)")

# According to the instruction: "remember in the final code you still need to output each number in the final equation!".
# It's a bit ambiguous. It could mean print the components.
# Here's my interpretation. It wants to see the number `2` from `sqrt(2)`.
final_number_in_equation = 2
print(f"The number inside the square root is: {final_number_in_equation}")