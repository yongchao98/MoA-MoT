from itertools import combinations

boxes = [18, 93, 23, 70, 72, 22, 96, 96]
lifters = [60, 46, 61]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    remaining_lifters = lifters[:]
    used_lifters = set()
    
    for box in boxes[:]:
        # Find the best combination of lifters to lift the box
        best_combination = None
        for r in range(1, len(remaining_lifters) + 1):
            for combo in combinations(enumerate(remaining_lifters), r):
                indices, capacities = zip(*combo)
                if sum(capacities) >= box and not used_lifters.intersection(indices):
                    if best_combination is None or sum(capacities) < sum(best_combination[1]):
                        best_combination = (indices, capacities)
        
        if best_combination:
            indices, _ = best_combination
            step.append((box, list(indices)))
            for index in indices:
                used_lifters.add(index)
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
for idx, step in enumerate(steps):
    print(f"Step {idx + 1}: {step}")