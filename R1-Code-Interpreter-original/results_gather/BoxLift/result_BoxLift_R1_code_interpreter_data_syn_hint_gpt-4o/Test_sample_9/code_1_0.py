from itertools import combinations

boxes = [355, 40, 260, 256, 250, 107, 301, 298, 72, 316, 50, 319, 243, 51, 393, 268, 141, 178, 270, 158, 354, 336, 55, 84]
lifters = [100, 90, 66, 141, 184, 78, 50]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return list(combo)
    return None

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = []
        remaining_weight = box
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifter
                used_lifters[j] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the box from the original list
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")