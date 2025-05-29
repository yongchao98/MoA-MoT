def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = lifters[:]
        used_lifters = [False] * len(lifters)
        
        for box in boxes[:]:
            for i in range(len(available_lifters)):
                if used_lifters[i]:
                    continue
                if available_lifters[i] >= box:
                    step.append((box, [i]))
                    used_lifters[i] = True
                    boxes.remove(box)
                    break
            else:
                # Try to combine lifters
                total_capacity = 0
                lifter_indices = []
                for i in range(len(available_lifters)):
                    if not used_lifters[i]:
                        total_capacity += available_lifters[i]
                        lifter_indices.append(i)
                        if total_capacity >= box:
                            step.append((box, lifter_indices))
                            for idx in lifter_indices:
                                used_lifters[idx] = True
                            boxes.remove(box)
                            break
        
        steps.append(step)
        if len(steps) > 7:
            return "Cannot lift all boxes in 7 steps or less."
    
    return steps

boxes = [190, 107, 357, 385, 343, 392, 97, 277, 231, 77, 147, 138, 371, 303, 173, 397, 123, 72, 241, 291]
lifters = [100, 70, 149, 198, 118, 100, 153]

steps = assign_lifters_to_boxes(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)