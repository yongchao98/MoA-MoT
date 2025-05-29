# Box weights and lifter capacities
boxes = [18, 27, 65, 31, 60, 73, 30, 67]
lifters = [66, 72, 43]

# Sort boxes in descending order to prioritize heavier boxes
boxes.sort(reverse=True)

# Sort lifters in descending order to prioritize stronger lifters
lifters.sort(reverse=True)

# Initialize steps
steps = []

# Try to lift all boxes in 4 steps or less
for step in range(4):
    step_lift = []
    available_lifters = lifters[:]
    
    for box in boxes[:]:
        lifter_indices = []
        total_capacity = 0
        
        # Try to find the smallest number of lifters to lift the box
        for i, capacity in enumerate(available_lifters):
            if total_capacity < box:
                lifter_indices.append(i)
                total_capacity += capacity
                available_lifters[i] = 0  # Mark lifter as used for this step
            
            if total_capacity >= box:
                break
        
        if total_capacity >= box:
            step_lift.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step_lift)
    
    if not boxes:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")