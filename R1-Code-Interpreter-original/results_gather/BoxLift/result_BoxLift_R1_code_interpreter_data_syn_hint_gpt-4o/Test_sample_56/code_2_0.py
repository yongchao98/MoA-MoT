from itertools import combinations

boxes = [21, 25, 138, 168, 105, 127, 154, 55, 155, 125, 121, 58]
lifters = [78, 98, 51, 75, 43]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Sort lifters in descending order
lifters.sort(reverse=True)

def find_lifter_combination(box, lifters, used_lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used_lifters[i] for i in combo):
                if sum(lifters[i] for i in combo) >= box:
                    for i in combo:
                        used_lifters[i] = True
                    return list(combo)
    return None

steps = []
used_boxes = [False] * len(boxes)

for step in range(5):
    current_step = []
    used_lifters = [False] * len(lifters)
    
    for i, box in enumerate(boxes):
        if used_boxes[i]:
            continue
        lifter_indices = find_lifter_combination(box, lifters, used_lifters)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            used_boxes[i] = True
    
    steps.append(current_step)
    if all(used_boxes):
        break

print(steps)