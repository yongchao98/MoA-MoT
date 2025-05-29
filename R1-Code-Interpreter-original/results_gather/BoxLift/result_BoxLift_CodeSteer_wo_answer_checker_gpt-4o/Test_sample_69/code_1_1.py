boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

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
                remaining_capacity[idx] = 0  # Mark lifter as used in this step
    
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")