boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
lifters = [184, 149, 64, 190, 59, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
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

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")