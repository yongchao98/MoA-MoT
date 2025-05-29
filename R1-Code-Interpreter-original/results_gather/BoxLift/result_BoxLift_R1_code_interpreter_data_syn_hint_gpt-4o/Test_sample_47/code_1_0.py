boxes = [280, 32, 168, 120, 112, 285, 46, 33, 148, 74, 295, 279, 207, 138, 153, 63]
lifters = [53, 92, 104, 86, 106, 68]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        # Try to find a combination of lifters to lift the box
        current_lifters = []
        remaining_capacity = box
        for i, lifter in enumerate(remaining_lifters):
            if not used_lifters[i] and lifter <= remaining_capacity:
                current_lifters.append(i)
                remaining_capacity -= lifter
                used_lifters[i] = True
            if remaining_capacity <= 0:
                break
        
        if remaining_capacity <= 0:
            # Box can be lifted
            step.append((box, current_lifters))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")