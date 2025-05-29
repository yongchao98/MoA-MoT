boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        total_capacity = 0
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j]:
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
output = ""
for i, step in enumerate(steps):
    if step:
        output += f"Step {i+1}: {step}\n"

print(f"<<<{output.strip()}>>>")