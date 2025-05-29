from itertools import combinations
from functools import lru_cache

# Box weights and lifter capacities
boxes = [355, 40, 260, 256, 250, 107, 301, 298, 72, 316, 50, 319, 243, 51, 393, 268, 141, 178, 270, 158, 354, 336, 55, 84]
lifters = [100, 90, 66, 141, 184, 78, 50]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

@lru_cache(None)
def backtrack(remaining_boxes, available_lifters, step):
    if not remaining_boxes:
        return step, []

    if step >= 9:
        return float('inf'), []

    min_steps = float('inf')
    best_assignment = []

    for box in remaining_boxes:
        for r in range(1, len(available_lifters) + 1):
            for combo in combinations(enumerate(available_lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= box:
                    new_lifters = tuple(l for i, l in enumerate(available_lifters) if i not in indices)
                    new_boxes = tuple(b for b in remaining_boxes if b != box)
                    steps, assignment = backtrack(new_boxes, new_lifters, step + 1)
                    if steps < min_steps:
                        min_steps = steps
                        best_assignment = [(box, list(indices))] + assignment
                    break
            else:
                continue
            break

    return min_steps, best_assignment

# Start backtracking
min_steps, assignment = backtrack(tuple(boxes), tuple(lifters), 0)

# Check if a valid solution was found
if min_steps == float('inf'):
    print("<<<No valid solution found within 9 steps>>>")
else:
    # Organize the assignment into steps
    steps = [[] for _ in range(min_steps)]
    for box, lifter_indices in assignment:
        for step in steps:
            if all(lifter not in [l for _, l in step] for lifter in lifter_indices):
                step.append((box, lifter_indices))
                break

    # Print the steps
    output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(steps))
    print(f"<<<{output}>>>")