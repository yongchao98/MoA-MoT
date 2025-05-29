from itertools import combinations

boxes = [196, 90, 142, 22, 143, 101, 43, 166, 79, 129, 77, 187]
lifters = [71, 58, 104, 105]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
while boxes and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    remaining_boxes = boxes[:]
    
    for box in remaining_boxes:
        # Try to find a combination of lifters to lift the box
        lifter_indices = []
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if not any(used_lifters[i] for i in combo):
                    if sum(lifters[i] for i in combo) >= box:
                        lifter_indices = list(combo)
                        for i in combo:
                            used_lifters[i] = True
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box, lifter_indices))
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")