boxes = [174, 117, 74, 120, 98, 65, 112, 61, 46, 87, 176, 110]
lifters = [119, 104, 67, 98]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 5:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        lifter_indices = []
        total_capacity = 0
        
        for i, capacity in enumerate(lifters):
            if not used_lifters[i] and total_capacity + capacity <= box:
                total_capacity += capacity
                lifter_indices.append(i)
                used_lifters[i] = True
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")