boxes = [34, 92, 289, 35, 140, 174, 244, 136, 80, 252, 136, 181, 92, 235, 61, 286]
lifters = [61, 60, 150, 96, 58, 91]

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

for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")