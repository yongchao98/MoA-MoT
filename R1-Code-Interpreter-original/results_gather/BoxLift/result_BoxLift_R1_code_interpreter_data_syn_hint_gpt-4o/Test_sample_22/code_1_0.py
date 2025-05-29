boxes = [215, 113, 62, 255, 246, 148, 283, 149, 295, 250, 114, 48, 245, 145, 153, 50]
lifters = [82, 55, 159, 56, 88]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    remaining_lifters.sort(reverse=True)
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        total_capacity = 0
        lifter_indices = []
        for i, lifter in enumerate(remaining_lifters):
            if not used_lifters[i] and total_capacity < box:
                lifter_indices.append(i)
                total_capacity += lifter
            if total_capacity >= box:
                break
        if total_capacity >= box:
            # Assign the box to the lifters
            step.append((box, lifter_indices))
            # Mark lifters as used
            for i in lifter_indices:
                used_lifters[i] = True
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")