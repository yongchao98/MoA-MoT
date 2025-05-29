boxes = [156, 123, 139, 191, 91, 193, 100, 47, 36, 169, 90, 117]
lifters = [52, 55, 63, 44]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to check if a combination of lifters can lift a box
def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, capacity in enumerate(lifters):
        if not used[i] and remaining_weight > 0:
            if capacity <= remaining_weight:
                lifter_indices.append(i)
                remaining_weight -= capacity
                used[i] = True
    return remaining_weight <= 0, lifter_indices

# Backtracking function to find the solution
def backtrack(boxes, lifters, step, steps, used):
    if not boxes:
        steps.append(step[:])
        return True
    if len(steps) >= 9:
        return False
    
    for i, box in enumerate(boxes):
        used_copy = used[:]
        can_lift_box, lifter_indices = can_lift(box, lifters, used_copy)
        if can_lift_box:
            step.append((box, lifter_indices))
            remaining_boxes = boxes[:i] + boxes[i+1:]
            if backtrack(remaining_boxes, lifters, step, steps, used_copy):
                return True
            step.pop()
    return False

# Initialize variables
steps = []
used = [False] * len(lifters)

# Start backtracking
backtrack(boxes, lifters, [], steps, used)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")