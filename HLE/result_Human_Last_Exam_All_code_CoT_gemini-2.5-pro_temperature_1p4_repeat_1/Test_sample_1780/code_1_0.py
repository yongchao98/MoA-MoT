import math

# Step 1: Identify the letters that can be freely arranged based on the connection rule.
# The letter 'W' cannot connect to any other letter, so it must be at the end.
# The letters 'L', 'N', 'S' can all connect to each other and to 'W'.
# Therefore, we need to find the number of ways to arrange 'L', 'N', and 'S'.
letters_to_arrange = ['L', 'N', 'S']

# Step 2: The number of ways is the number of permutations of these letters,
# which is n! (n-factorial).
n = len(letters_to_arrange)
result = math.factorial(n)

# Step 3: Create the equation string as requested.
# We want to show "3 * 2 * 1 = 6".
equation_parts = list(range(n, 0, -1))
equation_str = " * ".join(map(str, equation_parts))

# Step 4: Print the final result in the specified format.
print(f"The number of arrangements is the number of permutations of {n} letters.")
print(f"This is calculated as {n}! = {equation_str} = {result}.")
