from itertools import combinations

boxes = [256, 53, 286, 205, 261, 209, 242, 143, 354, 55, 72, 161, 268, 305, 43, 246, 146, 331, 84, 128, 88, 358, 71, 93]
lifters = [189, 58, 58, 79, 123, 192, 116]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    for box in boxes[:]:
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, selected_lifters = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(selected_lifters) >= box:
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    boxes.remove(box)
                    break
            if box not in boxes:
                break
    steps.append(step)
    if len(steps) > 8:
        break

# Print the steps
for i, step in enumerate(steps):
    print(f"Step {i+1}: {step}")