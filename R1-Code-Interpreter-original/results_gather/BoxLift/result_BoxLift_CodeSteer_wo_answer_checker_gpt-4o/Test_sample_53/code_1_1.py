boxes = [108, 206, 326, 343, 45, 374, 227, 366, 110, 136, 343, 298, 201, 168, 191, 243, 365, 375, 294, 317]
lifters = [147, 80, 144, 133, 69, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 9:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_capacity = lifters[:]
    
    for box in boxes[:]:
        lifter_indices = []
        total_capacity = 0
        
        for j, lifter in enumerate(remaining_capacity):
            if not used_lifters[j] and total_capacity < box:
                total_capacity += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
                remaining_capacity[j] = 0  # Mark this lifter as used in this step
            
            if total_capacity >= box:
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