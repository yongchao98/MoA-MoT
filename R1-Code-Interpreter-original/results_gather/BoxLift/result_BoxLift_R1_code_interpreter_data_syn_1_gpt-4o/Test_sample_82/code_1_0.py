boxes = [99, 159, 55, 86, 215, 56, 271, 167, 212, 87, 78, 211, 246, 153, 227, 287]
lifters = [143, 50, 70, 79, 136]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    remaining_lifters = lifters[:]
    used_lifters_indices = set()
    for box in boxes[:]:
        # Try to find a combination of lifters that can lift the box
        used_lifters = []
        total_capacity = 0
        for i, lifter in enumerate(remaining_lifters):
            if i not in used_lifters_indices and total_capacity < box:
                used_lifters.append(i)
                total_capacity += lifter
            if total_capacity >= box:
                break
        if total_capacity >= box:
            # Assign the lifters to this box
            step.append((box, used_lifters))
            # Mark the used lifters
            used_lifters_indices.update(used_lifters)
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")