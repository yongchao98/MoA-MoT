from itertools import combinations

# Box weights and lifter capacities
boxes = [49, 239, 263, 284, 45, 159, 306, 398, 279, 137, 218, 387, 161, 287, 171, 345, 373, 271, 60, 396]
lifters = [151, 91, 70, 108, 182, 200]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes:
    step = []
    used_lifters = [False] * len(lifters)
    
    for box in boxes[:]:
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, capacities = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(capacities) >= box:
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    boxes.remove(box)
                    break
            if box not in boxes:
                break
    
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")