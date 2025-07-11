import numpy as np
from scipy.integrate import quad

# Define the first simplified term of the integrand
def first_term(x):
    return 2**(-1/16) * (x**4) / np.sqrt(512 - x**8)

# Define the second simplified term of the integrand
def second_term(x):
    return 2**(1/16) * (x / np.sqrt(4 + x**2))**(1/4)

# Calculate the definite integral of the first term from 0 to 2
integral_1, error_1 = quad(first_term, 0, 2)

# Calculate the definite integral of the second term from 0 to 2
integral_2, error_2 = quad(second_term, 0, 2)

# The total integral is the sum of the two parts
total_integral = integral_1 + integral_2

# Output the results of the calculation
print(f"The integral of the first term is: {integral_1}")
print(f"The integral of the second term is: {integral_2}")
print(f"The value of the definite integral is the sum of the two parts.")
print(f"{integral_1} + {integral_2} = {total_integral}")
print(f"The final result is {round(total_integral)}.")
