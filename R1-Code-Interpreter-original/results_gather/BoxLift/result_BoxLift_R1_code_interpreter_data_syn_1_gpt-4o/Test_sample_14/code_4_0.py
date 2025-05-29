boxes = [364, 138, 172, 200, 44, 194, 195, 282, 189, 400, 233, 307, 71, 300, 306, 346, 61, 136, 168, 392, 97, 367, 118, 354]
lifters = [54, 124, 186, 103, 122, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 10:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Use a copy of the list to modify it during iteration
        lifter_indices = []
        remaining_weight = box
        
        # Try to find a combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j]:
                lifter_indices.append(j)
                remaining_weight -= lifter
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            # Mark the lifters as used only if the box can be lifted
            for index in lifter_indices:
                used_lifters[index] = True
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the box once it is lifted
    
    steps.append(step)

print(steps)