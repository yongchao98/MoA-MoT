from itertools import combinations

boxes = [14, 13, 59, 47, 73, 64, 58, 75]
lifters = [74, 57, 51]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

for step in range(4):
    current_step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        
        # Try to find a single lifter or combination of lifters to lift the box
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and lifter >= box:
                current_step.append((box, [j]))
                used_boxes[i] = True
                used_lifters[j] = True
                break
        else:
            # Try to find a combination of lifters
            for comb in combinations(range(len(lifters)), 2):
                if not used_lifters[comb[0]] and not used_lifters[comb[1]]:
                    if lifters[comb[0]] + lifters[comb[1]] >= box:
                        current_step.append((box, list(comb)))
                        used_boxes[i] = True
                        used_lifters[comb[0]] = True
                        used_lifters[comb[1]] = True
                        break
    
    steps.append(current_step)
    if all(used_boxes):
        break

print(steps)