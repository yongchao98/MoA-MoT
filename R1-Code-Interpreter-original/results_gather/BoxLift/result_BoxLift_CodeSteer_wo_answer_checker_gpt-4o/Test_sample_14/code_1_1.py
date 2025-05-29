from itertools import combinations

boxes = [364, 138, 172, 200, 44, 194, 195, 282, 189, 400, 233, 307, 71, 300, 306, 346, 61, 136, 168, 392, 97, 367, 118, 354]
lifters = [54, 124, 186, 103, 122, 112]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 10:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for box in boxes:
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(enumerate(lifters), r):
                indices, capacities = zip(*combo)
                if not any(used_lifters[i] for i in indices) and sum(capacities) >= box:
                    step.append((box, list(indices)))
                    for i in indices:
                        used_lifters[i] = True
                    remaining_boxes.remove(box)
                    break
            else:
                continue
            break
    
    boxes = remaining_boxes
    steps.append(step)

# Format the output
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")