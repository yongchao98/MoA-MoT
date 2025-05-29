from itertools import combinations

boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
lifters = [200, 91, 128, 198, 196, 119, 188]

# Sort boxes in descending order
boxes.sort(reverse=True)

steps = []
box_index = 0

while box_index < len(boxes) and len(steps) < 6:
    step = []
    used_lifters = [False] * len(lifters)
    
    for i in range(len(boxes)):
        if box_index >= len(boxes):
            break
        box_weight = boxes[box_index]
        lifter_indices = []
        
        # Try to find a combination of lifters that can lift the box
        for r in range(1, len(lifters) + 1):
            for combo in combinations(range(len(lifters)), r):
                if all(not used_lifters[j] for j in combo):
                    total_capacity = sum(lifters[j] for j in combo)
                    if total_capacity >= box_weight:
                        lifter_indices = list(combo)
                        for j in combo:
                            used_lifters[j] = True
                        break
            if lifter_indices:
                break
        
        if lifter_indices:
            step.append((box_weight, lifter_indices))
            box_index += 1
    
    steps.append(step)

print(steps)