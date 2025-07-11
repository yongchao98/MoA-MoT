import sys

# Define the given parameters
c = 0.95  # Consistency level of identifier resolution
b = 3     # Branching factor of semantic version control

# --- Define the FAIR score components based on the plan ---

# Findability (f) is limited by the identifier resolution consistency.
f_max = c

# Accessibility (a) is also limited by the identifier resolution consistency,
# as data cannot be accessed if its identifier cannot be resolved.
a_max = c

# Interoperability (i) is inversely affected by the branching factor 'b'.
# A higher branching factor creates more versions to manage, reducing interoperability.
# The chosen model is i = 1/b.
i_max = 1 / b

# Reusability (r) is also inversely affected by 'b' for the same reasons.
# The chosen model is r = 1/b.
r_max = 1 / b

# --- Calculate the theoretical maximum FAIR score R ---

# The overall FAIR score R is the average of its four components.
R_max = (f_max + a_max + i_max + r_max) / 4

# --- Output the result ---

# We need to explicitly display all the numbers in the final equation.
# Let's format the floating point numbers for readability.
f_str = f"{f_max:.2f}"
a_str = f"{a_max:.2f}"
i_str = f"{i_max:.4f}"
r_str = f"{r_max:.4f}"
R_str = f"{R_max:.4f}"

# The final print statement constructs the full equation as a string.
print("Based on the provided parameters and a quantitative FAIR model:")
print(f"Findability (f) = {f_str}")
print(f"Accessibility (a) = {a_str}")
print(f"Interoperability (i) = 1/{b} = {i_str}")
print(f"Reusability (r) = 1/{b} = {r_str}")
print("\nThe theoretical maximum value of R is calculated as:")
print(f"R = (f + a + i + r) / 4")
print(f"R = ({f_str} + {a_str} + {i_str} + {r_str}) / 4")
# To show the sum before the division for clarity:
numerator_sum = f_max + a_max + i_max + r_max
print(f"R = {numerator_sum:.4f} / 4")
print(f"R = {R_str}")

# Add the final answer tag for parsing. The value 77/120 = 0.641666... is rounded to four decimal places.
sys.stdout.write("<<<0.6417>>>")