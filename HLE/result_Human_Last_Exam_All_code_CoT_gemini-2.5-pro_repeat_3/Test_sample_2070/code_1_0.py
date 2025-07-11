import sys

# Check if the Python version is 3.13.1, as specified.
# This part is for ensuring compatibility but is not part of the core logic.
if sys.version_info.major != 3 or sys.version_info.minor != 13 or sys.version_info.micro != 1:
    print("Warning: This code was designed for Python 3.13.1. Output may vary on other versions.", file=sys.stderr)

# The logic as determined in the step-by-step analysis.
m = True
a = set((lambda f: [f.append(f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m, m]))
b = set((lambda f: [f.append(f[~(m ^ -~m) ^ m] + f[-m] + f[--~m]) or f[---m] for _ in range(len("1 2 3 4 5"))] and f)([not m]*(m << m)+[m]))

# Perform the final operation
result_set = b & a
final_value = sum(result_set)

# Get the numbers from the resulting set to display the equation
sorted_numbers = sorted(list(result_set))
equation_str = " + ".join(map(str, sorted_numbers))

# Print the final equation as requested
print(f"{equation_str} = {final_value}")
