from itertools import combinations

boxes = [50, 219, 120, 245, 204, 172, 274, 225, 94, 270, 158, 253, 174, 40, 252, 124]
lifters = [54, 61, 106, 60, 63]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        total_capacity = 0
        for i, lifter in enumerate(remaining_lifters):
            if not used_lifters[i] and total_capacity < box:
                lifter_indices.append(i)
                total_capacity += lifter
                used_lifters[i] = True
            if total_capacity >= box:
                break
        if total_capacity >= box:
            # Assign the box to the lifters
            step.append((box, lifter_indices))
            # Remove the box from the list
            boxes.remove(box)
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")