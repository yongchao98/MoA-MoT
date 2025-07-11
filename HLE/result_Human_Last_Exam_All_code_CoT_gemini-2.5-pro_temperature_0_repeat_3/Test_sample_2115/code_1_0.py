import math

# The problem reduces to calculating the definite integral of the initial condition
# on the slice z = -y, from x=0 to x=1.
# The analytical solution to the integral is 3 * (ln(3) - ln(e^2 + e + 1)).

# Calculate the components of the final expression.
e = math.exp(1)
val_e_squared = e**2
term_in_log = val_e_squared + e + 1
log_of_3 = math.log(3)
log_of_term = math.log(term_in_log)

# Calculate the final result
final_result = 3 * (log_of_3 - log_of_term)

# Print the breakdown of the calculation as requested.
print("The integral is evaluated as: 3 * (ln(3) - ln(e^2 + e + 1))")
print("\n--- Numerical Breakdown ---")
print(f"e = {e}")
print(f"e^2 + e + 1 = {term_in_log}")
print(f"ln(3) = {log_of_3}")
print(f"ln(e^2 + e + 1) = {log_of_term}")
print("\n--- Final Equation ---")
print(f"Result = 3 * ({log_of_3} - {log_of_term})")
print(f"Result = {final_result}")
