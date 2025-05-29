from itertools import combinations

boxes = [299, 338, 270, 207, 76, 195, 101, 89, 277, 236, 240, 342, 47, 207, 220, 349, 185, 124, 309, 112, 266, 278, 162, 119]
lifters = [198, 55, 76, 137, 163, 117, 146]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Precompute feasible combinations of lifters for each box
feasible_combinations = {}
for box in boxes:
    feasible_combinations[box] = []
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, weights = zip(*combo)
            if sum(weights) >= box:
                feasible_combinations[box].append(indices)

# Backtracking to find the solution
def backtrack(remaining_boxes, used_lifters, steps):
    if not remaining_boxes:
        return steps

    box = remaining_boxes[0]
    for combo in feasible_combinations[box]:
        if all(not used_lifters[i] for i in combo):
            # Mark lifters as used
            for i in combo:
                used_lifters[i] = True

            # Recurse with the remaining boxes
            result = backtrack(remaining_boxes[1:], used_lifters, steps + [(box, list(combo))])
            if result:
                return result

            # Unmark lifters
            for i in combo:
                used_lifters[i] = False

    return None

# Initialize used lifters and find the solution
used_lifters = [False] * len(lifters)
solution = backtrack(boxes, used_lifters, [])

# Output the steps
if solution:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within the constraints.")