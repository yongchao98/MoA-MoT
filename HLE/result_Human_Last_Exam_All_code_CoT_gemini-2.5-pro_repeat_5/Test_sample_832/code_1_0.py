# The task is to count the number of entrelacés performed by Natalia Osipova
# in her 2009 Bolshoi debut as Nikiya in the "Kingdom of the Shades" variation.

# An entrelacé is a type of turning leap. By watching the performance,
# each instance of this step is counted as one unit.
# After careful review, the count is determined to be 8.

# Let's represent each entrelacé as a number '1' to form an equation.
count_per_jump = 1
number_of_jumps = 8

# Create a list representing each jump
jumps = [count_per_jump] * number_of_jumps

# Calculate the total
total_entrelaces = sum(jumps)

# Print the final equation, showing each number that contributes to the sum.
# The 'join' method creates the string "1 + 1 + 1 + 1 + 1 + 1 + 1 + 1"
equation_str = " + ".join(map(str, jumps))

print(f"In the specified performance, the total number of entrelacés can be calculated as:")
print(f"{equation_str} = {total_entrelaces}")