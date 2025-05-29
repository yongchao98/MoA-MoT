boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
lifters = [144, 152, 155, 109, 115, 152, 162]

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
            if not used_lifters[j]:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity >= box:
                    step.append((box, lifter_indices))
                    for idx in lifter_indices:
                        used_lifters[idx] = True
                    boxes.remove(box)
                    break
        
        # If the box cannot be lifted by any combination of lifters, raise an error
        if total_capacity < box:
            raise ValueError(f"Box with weight {box} cannot be lifted by the available lifters.")
    
    steps.append(step)
    if len(steps) > 6:
        raise ValueError("Cannot lift all boxes in 6 steps or less.")

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output.strip()}>>>")