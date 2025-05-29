def assign_lifters_to_boxes(boxes, lifters):
    boxes = sorted(boxes, reverse=True)
    lifters = sorted(lifters, reverse=True)
    steps = []

    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for i, box in enumerate(boxes):
            if not boxes:
                break
            lifter_indices = []
            remaining_weight = box
            
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and lifter <= remaining_weight:
                    lifter_indices.append(j)
                    remaining_weight -= lifter
                    used_lifters[j] = True
                if remaining_weight <= 0:
                    break
            
            if remaining_weight <= 0:
                step.append((box, lifter_indices))
                boxes[i] = None
        
        boxes = [box for box in boxes if box is not None]
        steps.append(step)
    
    return steps

boxes = [351, 64, 363, 78, 283, 68, 211, 200, 288, 295, 319, 271, 305, 219, 310, 365, 390, 56, 396, 332, 130, 136, 221, 386]
lifters = [134, 184, 163, 144, 109, 179, 64]

steps = assign_lifters_to_boxes(boxes, lifters)

output = ""
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"

print(f"<<<{output}>>>")