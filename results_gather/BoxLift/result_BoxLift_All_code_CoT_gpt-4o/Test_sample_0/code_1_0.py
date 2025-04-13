boxes = [43, 88, 84, 64, 56, 87, 20, 57, 82, 14]
lifters = [51, 49, 48]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and remaining_weight > 0:
                if lifter >= remaining_weight:
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    remaining_weight = 0
                    break
                elif lifter < remaining_weight:
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    remaining_weight -= lifter
        
        if remaining_weight == 0:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")