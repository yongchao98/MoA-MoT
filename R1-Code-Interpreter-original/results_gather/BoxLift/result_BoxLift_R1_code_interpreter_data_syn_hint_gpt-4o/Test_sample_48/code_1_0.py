from itertools import combinations

boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
lifters = [144, 152, 155, 109, 115, 152, 162]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
remaining_boxes = boxes[:]

while remaining_boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in remaining_boxes[:]:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        total_capacity = 0
        
        # Check all combinations of lifters
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo):
                    total_capacity = sum(lifters[i] for i in combo)
                    if total_capacity >= box:
                        lifter_indices = list(combo)
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            for i in lifter_indices:
                used_lifters[i] = True
            remaining_boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")