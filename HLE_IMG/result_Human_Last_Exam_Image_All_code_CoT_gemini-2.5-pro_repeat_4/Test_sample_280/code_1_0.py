# The program's visual structure is based on a repeating pattern of colored bars.
# The number of distinct colors in the repeating pattern.
colors_in_pattern = 6

# In the common version of this puzzle, the pattern is repeated 7 times.
repetitions = 7

# The total number of columns, which is the likely intended answer.
total_columns = colors_in_pattern * repetitions

# The problem asks to output the final equation.
print(f"{colors_in_pattern} * {repetitions} = {total_columns}")