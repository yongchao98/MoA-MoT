boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_capacity = lifters[:]
    
    for box in boxes[:]:
        lifter_indices = []
        total_capacity = 0
        
        for i, capacity in enumerate(remaining_capacity):
            if not used_lifters[i] and total_capacity + capacity <= box:
                total_capacity += capacity
                lifter_indices.append(i)
                used_lifters[i] = True
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
            for idx in lifter_indices:
                remaining_capacity[idx] = 0  # Mark these lifters as used for this step
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")