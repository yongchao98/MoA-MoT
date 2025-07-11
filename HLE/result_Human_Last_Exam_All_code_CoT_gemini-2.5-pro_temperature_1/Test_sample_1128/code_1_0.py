import math

# The number of symbols, n, is not specified in the problem.
# We will use n=5 as an example to demonstrate the calculation.
n = 5

# The number of cool strings of maximal length is n!
# We check if n is a non-negative integer first.
if not isinstance(n, int) or n < 0:
  print("Error: n must be a non-negative integer.")
else:
  # Calculate n factorial
  result = math.factorial(n)

  # Explain the formula
  print(f"For n = {n} symbols, the number of cool strings of maximal length ({3*n}) is n!.")
  
  # Create the string for the full equation
  if n == 0:
    equation_str = "1"
  else:
    # Generate parts of the equation: e.g., ["5", "4", "3", "2", "1"]
    equation_parts = [str(i) for i in range(n, 0, -1)]
    equation_str = " * ".join(equation_parts)
  
  # Print the final formatted output
  print(f"The calculation is: {n}! = {equation_str} = {result}")
