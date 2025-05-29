boxes = [351, 64, 363, 78, 283, 68, 211, 200, 288, 295, 319, 271, 305, 219, 310, 365, 390, 56, 396, 332, 130, 136, 221, 386]
lifters = [134, 184, 163, 144, 109, 179, 64]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        total_capacity = 0
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity < box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
            
            if total_capacity >= box:
                break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")