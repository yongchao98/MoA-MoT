from itertools import combinations

boxes = [159, 197, 211, 302, 51, 251, 299, 229, 384, 196, 377, 283, 127, 312, 228, 280, 153, 199, 169, 151]
lifters = [52, 80, 175, 75, 136, 143]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        for r in range(1, len(lifters) + 1):
            for lifter_indices in combinations(range(len(lifters)), r):
                if all(not used_lifters[i] for i in lifter_indices):
                    total_capacity = sum(lifters[i] for i in lifter_indices)
                    if total_capacity >= box:
                        step.append((box, list(lifter_indices)))
                        for i in lifter_indices:
                            used_lifters[i] = True
                        boxes.remove(box)
                        break
            if box not in boxes:
                break
    
    steps.append(step)
    if len(steps) > 10:
        print("Unable to lift all boxes in 10 steps or less.")
        break

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")