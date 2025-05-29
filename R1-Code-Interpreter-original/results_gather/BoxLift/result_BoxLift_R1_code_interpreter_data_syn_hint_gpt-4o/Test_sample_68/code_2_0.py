from itertools import combinations

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
        
        # Try all combinations of lifters
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= remaining_weight and all(not used_lifters[i] for i in indices):
                    lifter_indices = list(indices)
                    for i in indices:
                        used_lifters[i] = True
                    break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")