boxes = [190, 107, 357, 385, 343, 392, 97, 277, 231, 77, 147, 138, 371, 303, 173, 397, 123, 72, 241, 291]
lifters = [100, 70, 149, 198, 118, 100, 153]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes and len(steps) < 7:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = []
        remaining_weight = box
        
        # Sort lifters by capacity in descending order
        lifter_pairs = sorted(enumerate(lifters), key=lambda x: x[1], reverse=True)
        
        for j, lifter in lifter_pairs:
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifter
                used_lifters[j] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")