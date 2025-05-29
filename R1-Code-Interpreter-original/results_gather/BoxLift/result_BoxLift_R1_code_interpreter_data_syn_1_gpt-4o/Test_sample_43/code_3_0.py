from itertools import combinations

boxes = [247, 236, 255, 64, 58, 144, 218, 200, 320, 87, 256, 111, 379, 318, 367, 88, 136, 125, 172, 93]
lifters = [63, 51, 130, 110, 117, 132]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:  # Iterate over a copy of the list
        lifter_indices = []
        remaining_weight = box
        # Try all combinations of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in combo):
                    total_capacity = sum(lifters[i] for i in combo)
                    if total_capacity >= remaining_weight:
                        lifter_indices = list(combo)
                        for i in combo:
                            used_lifters[i] = True
                        break
            if lifter_indices:
                break
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)  # Remove the box once it's lifted
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i + 1}: {step}")