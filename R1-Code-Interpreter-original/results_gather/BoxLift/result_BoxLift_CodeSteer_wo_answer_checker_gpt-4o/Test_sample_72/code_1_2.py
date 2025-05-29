from itertools import combinations

boxes = [280, 109, 195, 129, 54, 216, 157, 193, 44, 233, 171, 218, 204, 120, 131, 151]
lifters = [49, 131, 66, 124, 91]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(enumerate(available_lifters), r):
            indices, capacities = zip(*combo)
            if sum(capacities) >= box:
                return list(indices)
    return None

def backtrack(remaining_boxes, available_lifters, step_count, steps):
    if not remaining_boxes:
        return step_count, steps

    if step_count > 8:
        return float('inf'), []

    state = (tuple(remaining_boxes), tuple(available_lifters))
    if state in memo:
        return memo[state]

    min_steps = float('inf')
    best_steps = []

    for i, box in enumerate(remaining_boxes):
        lifter_indices = can_lift(box, available_lifters)
        if lifter_indices is not None:
            new_lifters = [l for j, l in enumerate(available_lifters) if j not in lifter_indices]
            new_boxes = remaining_boxes[:i] + remaining_boxes[i+1:]
            current_step = steps + [(box, lifter_indices)]
            steps_taken, step_details = backtrack(new_boxes, new_lifters, step_count + 1, current_step)
            if steps_taken < min_steps:
                min_steps = steps_taken
                best_steps = step_details

    memo[state] = (min_steps, best_steps)
    return min_steps, best_steps

# Start backtracking
min_steps, detailed_steps = backtrack(boxes, lifters, 0, [])

# Output the result
if min_steps <= 8:
    output = []
    for i, step in enumerate(detailed_steps):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No valid solution within 8 steps.")