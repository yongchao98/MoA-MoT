from itertools import combinations

boxes = [156, 123, 139, 191, 91, 193, 100, 47, 36, 169, 90, 117]
lifters = [52, 55, 63, 44]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to check if a combination of lifters can lift a box
def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box:
                    return True, list(combo)
    return False, []

# Backtracking function to find the solution
def backtrack(boxes, lifters, step, steps, used, max_steps):
    if not boxes:
        steps.append(step[:])
        return True
    if len(steps) >= max_steps:
        return False
    
    for i, box in enumerate(boxes):
        can_lift_box, lifter_indices = can_lift(box, lifters, used)
        if can_lift_box:
            # Mark lifters as used
            for idx in lifter_indices:
                used[idx] = True
            step.append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step, steps, used, max_steps):
                return True
            # Backtrack: unmark lifters and remove the step
            for idx in lifter_indices:
                used[idx] = False
            step.pop()
    
    return False

# Initialize variables
steps = []
used = [False] * len(lifters)

# Start backtracking with a maximum of 9 steps
backtrack(boxes, lifters, [], steps, used, 9)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")