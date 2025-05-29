from itertools import combinations

# Box weights and lifter capacities
boxes = [187, 170, 257, 238, 229, 202, 239, 125, 197, 243, 110, 50, 48, 168, 42, 292]
lifters = [115, 55, 142, 79, 63]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a solution using backtracking
def backtrack(steps, current_step, remaining_boxes, available_lifters):
    if not remaining_boxes:
        return steps

    if current_step >= 10:
        return None

    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            new_step = []
            used_lifters = set()
            for box in remaining_boxes:
                if sum(combo) >= box:
                    new_step.append((box, [lifters.index(l) for l in combo]))
                    used_lifters.update(lifters.index(l) for l in combo)
                    remaining_boxes.remove(box)
                    break
            if new_step:
                new_steps = steps + [new_step]
                result = backtrack(new_steps, current_step + 1, remaining_boxes, [l for l in available_lifters if lifters.index(l) not in used_lifters])
                if result:
                    return result
    return None

# Start backtracking
solution = backtrack([], 0, boxes[:], lifters)

# Print the solution
if solution:
    output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(solution))
    print(f"<<<{output}>>>")
else:
    print("No solution found within 10 steps.")