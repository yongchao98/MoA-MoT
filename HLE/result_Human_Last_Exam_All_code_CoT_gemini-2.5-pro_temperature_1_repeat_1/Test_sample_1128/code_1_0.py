import math

# Step 1: Determine the value of n.
# The problem is stated for a general 'n', but the structure of the problem,
# particularly the maximal length being 3n, suggests using n=3.
n = 3

# Step 2: The number of cool strings of maximal length is n!
# We will calculate the factorial of n.

if n < 0:
    print("Factorial is not defined for negative numbers.")
elif n == 0:
    # The equation for 0! is simply 1.
    result = 1
    print(f"{n}! = {result}")
else:
    # Calculate the factorial and prepare the string for the equation.
    result = math.factorial(n)
    
    # Create a list of numbers from 1 to n for the equation string.
    factors = [str(i) for i in range(1, n + 1)]
    
    # Step 3: Print the final equation and the result.
    # The instruction is to output each number in the final equation.
    print(f"For n = {n}, the number of cool strings is n!")
    print(f"The calculation is:")
    print(f"{n}! = {' * '.join(factors)} = {result}")
