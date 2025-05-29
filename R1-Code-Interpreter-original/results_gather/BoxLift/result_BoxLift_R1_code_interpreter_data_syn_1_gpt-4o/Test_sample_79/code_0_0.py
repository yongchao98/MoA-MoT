boxes = [368, 220, 131, 287, 153, 175, 153, 217, 262, 115, 365, 394, 291, 222, 142, 126, 189, 82, 159, 367, 356, 168, 188, 355]
lifters = [200, 91, 128, 198, 196, 119, 188]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

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
        remaining_weight = box_weight
        
        for j in range(len(lifters)):
            if not used_lifters[j] and lifters[j] <= remaining_weight:
                lifter_indices.append(j)
                remaining_weight -= lifters[j]
                used_lifters[j] = True
                if remaining_weight <= 0:
                    break
        
        if remaining_weight <= 0:
            step.append((box_weight, lifter_indices))
            box_index += 1
    
    steps.append(step)

print(steps)