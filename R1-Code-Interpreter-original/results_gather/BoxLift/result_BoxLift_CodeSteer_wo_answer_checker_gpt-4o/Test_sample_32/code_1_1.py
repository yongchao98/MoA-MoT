def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []

    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for box in boxes[:]:
            total_capacity = 0
            lifter_indices = []
            
            for j, lifter in enumerate(lifters):
                if not used_lifters[j] and total_capacity + lifter <= box:
                    total_capacity += lifter
                    lifter_indices.append(j)
                    used_lifters[j] = True
                    if total_capacity == box:
                        break
            
            if total_capacity >= box:
                step.append((box, lifter_indices))
                boxes.remove(box)
        
        steps.append(step)
        if len(steps) > 11:
            break
    
    return steps

boxes = [257, 47, 186, 194, 237, 135, 72, 204, 268, 83, 242, 41, 214, 157, 115, 242]
lifters = [77, 41, 40, 54, 132]

steps = assign_lifters_to_boxes(boxes, lifters)

output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output}>>>")