import math

# Define the boundaries of the integral
x_upper = 1
x_lower = 0

# Calculate the term inside the logarithm at the upper boundary x=1
val_upper = math.exp(2 * x_upper) + math.exp(x_upper) + 1

# Calculate the term inside the logarithm at the lower boundary x=0
val_lower = math.exp(2 * x_lower) + math.exp(x_lower) + 1

# Calculate the natural logarithms of these values
log_val_upper = math.log(val_upper)
log_val_lower = math.log(val_lower)

# Calculate the final result of the definite integral
final_result = -3 * (log_val_upper - log_val_lower)

# Print the step-by-step evaluation of the final expression
print("The analytical solution to the integral is: -3 * [ln(e^2 + e + 1) - ln(e^0 + e^0 + 1)]")
print(f"This evaluates to: -3 * (ln({val_upper}) - ln({val_lower}))")
print(f"Using the numerical values for the logarithms, the final equation is:")
print(f"-3 * ({log_val_upper} - {log_val_lower}) = {final_result}")