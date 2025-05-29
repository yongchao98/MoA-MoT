boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
lifters = [71, 58, 104, 105]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for box in remaining_boxes:
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_weight -= lifter
            if remaining_weight <= 0:
                break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")