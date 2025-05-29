boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
lifters = [55, 85, 43, 100, 47]

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
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")