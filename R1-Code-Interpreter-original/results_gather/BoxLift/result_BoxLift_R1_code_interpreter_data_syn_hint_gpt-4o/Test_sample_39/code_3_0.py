from itertools import combinations

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
        # Try to find a single lifter or a combination of lifters to lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, combo_lifters = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(combo_lifters) >= box:
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    boxes.remove(box)
                    break
            else:
                continue
            break
    
    steps.append(step)

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")