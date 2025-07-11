import math

# This script calculates the expected maximum earthquake magnitude based on the derived formula.
# The formula is E = pi / (2 * log(2)) - 1.

# --- Components of the equation ---
pi_val = math.pi
log_2 = math.log(2)
two_log_2 = 2 * log_2
term1 = pi_val / two_log_2
result = term1 - 1

# --- Printing the step-by-step calculation ---
print("The final formula for the expected maximum magnitude is: E = pi / (2 * log(2)) - 1")
print(f"Value of pi: {pi_val}")
print(f"Value of log(2): {log_2}")
print(f"Value of 2 * log(2): {two_log_2}")
print(f"Value of pi / (2 * log(2)): {term1}")
print(f"Final Result (E): {result}")
