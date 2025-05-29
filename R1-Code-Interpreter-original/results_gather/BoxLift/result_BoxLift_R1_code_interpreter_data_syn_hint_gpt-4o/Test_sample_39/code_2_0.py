boxes = [256, 53, 286, 205, 261, 209, 242, 143, 354, 55, 72, 161, 268, 305, 43, 246, 146, 331, 84, 128, 88, 358, 71, 93]
lifters = [189, 58, 58, 79, 123, 192, 116]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 8:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        current_weight = 0
        lifter_indices = []
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and current_weight + lifter <= box:
                current_weight += lifter
                lifter_indices.append(j)
                used_lifters[j] = True
            if current_weight >= box:
                step.append((box, lifter_indices))
                boxes.remove(box)
                break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")