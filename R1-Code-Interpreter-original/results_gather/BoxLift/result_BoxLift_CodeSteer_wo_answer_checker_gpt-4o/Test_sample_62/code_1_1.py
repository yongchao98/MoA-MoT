boxes = [91, 207, 152, 47, 209, 54, 251, 176, 194, 221, 152, 141, 128, 159, 57, 184]
lifters = [149, 131, 113, 109, 124]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to find a valid assignment of lifters to boxes
def find_steps(boxes, lifters, step, steps, max_steps):
    if not boxes:
        return True
    if step >= max_steps:
        return False

    for i, box in enumerate(boxes):
        # Try to find a combination of lifters to lift the box
        for j in range(1 << len(lifters)):
            total_capacity = 0
            lifter_indices = []
            for k in range(len(lifters)):
                if j & (1 << k):
                    total_capacity += lifters[k]
                    lifter_indices.append(k)
            if total_capacity >= box:
                # Assign this combination to the current step
                steps[step].append((box, lifter_indices))
                # Remove the used lifters temporarily
                new_lifters = [lifters[k] for k in range(len(lifters)) if k not in lifter_indices]
                # Recurse to the next box
                if find_steps(boxes[:i] + boxes[i+1:], new_lifters, step + 1, steps, max_steps):
                    return True
                # Backtrack
                steps[step].pop()
    return False

# Initialize steps
steps = [[] for _ in range(6)]
if find_steps(boxes, lifters, 0, steps, 6):
    output = []
    for i, step in enumerate(steps):
        if step:
            output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")