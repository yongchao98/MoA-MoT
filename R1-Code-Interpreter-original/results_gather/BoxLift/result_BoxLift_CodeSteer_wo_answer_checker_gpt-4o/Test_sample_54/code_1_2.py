from itertools import combinations

boxes = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
lifters = [167, 156, 74, 177, 77, 60]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Precompute all possible combinations of lifters that can lift each box
lifter_combinations = []
for box in boxes:
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box:
                valid_combinations.append(indices)
    lifter_combinations.append(valid_combinations)

# Use a branch-and-bound strategy to find the minimum number of steps
def find_min_steps(box_index, used_lifters, steps):
    if box_index == len(boxes):
        return steps
    
    min_steps = float('inf')
    for combo in lifter_combinations[box_index]:
        if all(not used_lifters[i] for i in combo):
            for i in combo:
                used_lifters[i] = True
            min_steps = min(min_steps, find_min_steps(box_index + 1, used_lifters, steps + 1))
            for i in combo:
                used_lifters[i] = False
    
    return min_steps

# Initialize used lifters array
used_lifters = [False] * len(lifters)

# Find the minimum number of steps
min_steps = find_min_steps(0, used_lifters, 0)

# Output the result
print(f"Minimum number of steps: {min_steps}")