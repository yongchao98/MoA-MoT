# 1. Define the numbers derived from the haiku.
twice_fifteen = 30
august = 8
result = 1

# 2. Solve for the unknown variable 'x'.
# The riddle implies the equation: (30 - 8) / x = 1
# This simplifies to 22 / x = 1, so x must be 22.
x = twice_fifteen - august

# 3. Print the full equation as requested, showing each number.
print(f"The equation from the haiku is: ({twice_fifteen} - {august}) / {x} = {result}")

# 4. Convert the numerical answer to its corresponding letter of the alphabet.
# The 1st letter is 'A', 2nd is 'B', ..., 22nd is 'V'.
# We can find this by using the ASCII character code. 'A' starts at 65.
alphabetical_answer = chr(x + 64)

# 5. Print the final result.
print(f"The solution to the equation is {x}, which corresponds to the letter '{alphabetical_answer}'.")
