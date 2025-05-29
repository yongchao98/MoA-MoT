boxes = [299, 338, 270, 207, 76, 195, 101, 89, 277, 236, 240, 342, 47, 207, 220, 349, 185, 124, 309, 112, 266, 278, 162, 119]
lifters = [198, 55, 76, 137, 163, 117, 146]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        total_capacity = 0
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                if total_capacity == box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Output the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")