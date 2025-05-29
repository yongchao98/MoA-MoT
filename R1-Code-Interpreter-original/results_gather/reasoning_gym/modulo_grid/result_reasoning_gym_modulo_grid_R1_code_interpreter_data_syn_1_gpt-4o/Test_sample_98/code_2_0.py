# Known differences
differences = [13, 5, 13, 8, 8, 2, -5, 12, -2, 2, 5, -2, -8, -8, -13, 4, -5, -13, -5, -12]

# Function to find repeating cycles
def find_repeating_cycle(differences):
    n = len(differences)
    for cycle_length in range(1, n // 2 + 1):
        cycle = differences[:cycle_length]
        if cycle * (n // cycle_length) == differences[:cycle_length * (n // cycle_length)]:
            return cycle
    return None

# Find the repeating cycle
cycle = find_repeating_cycle(differences)

# Print the repeating cycle
print("Repeating Cycle:", cycle)