boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifters = [61, 45, 66, 38, 52]

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

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")