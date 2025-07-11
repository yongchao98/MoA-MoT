import math

n = 9
# The cardinality is given by n!
# We will calculate 9!
result = math.factorial(n)

# We want to print the final equation.
# First, let's represent the factorial calculation as a string
equation_str = "9! = " + " * ".join(map(str, range(1, n + 1)))

print(f"The cardinality of [Gamma(F)](*,*) for n={n} is n!.")
print(f"For n = {n}, we need to compute {n}!:")
print(equation_str + f" = {result}")
