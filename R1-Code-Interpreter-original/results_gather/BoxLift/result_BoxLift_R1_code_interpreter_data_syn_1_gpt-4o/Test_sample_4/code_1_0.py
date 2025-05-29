from itertools import combinations

boxes = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
lifters = [184, 149, 64, 190, 59, 51]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes[:]):  # Use a copy of the list to iterate
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        remaining_weight = box
        
        # Check all combinations of lifters
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo):
                    total_capacity = sum(lifters[j] for j in combo)
                    if total_capacity >= remaining_weight:
                        lifter_indices = list(combo)
                        for j in combo:
                            used_lifters[j] = True
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")