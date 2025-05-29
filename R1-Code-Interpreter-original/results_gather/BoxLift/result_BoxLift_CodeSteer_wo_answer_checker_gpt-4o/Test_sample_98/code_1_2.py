def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []
    max_steps = 6

    for _ in range(max_steps):
        step = []
        used_lifters = [False] * len(lifters)
        
        for box in boxes[:]:  # Iterate over a copy of the list
            lifter_indices = []
            remaining_weight = box
            
            for i, capacity in enumerate(lifters):
                if not used_lifters[i] and capacity >= remaining_weight:
                    lifter_indices.append(i)
                    used_lifters[i] = True
                    remaining_weight -= capacity
                    break
            
            if remaining_weight > 0:  # If a single lifter can't lift the box
                for i, capacity in enumerate(lifters):
                    if not used_lifters[i] and capacity <= remaining_weight:
                        lifter_indices.append(i)
                        used_lifters[i] = True
                        remaining_weight -= capacity
                        if remaining_weight <= 0:
                            break
            
            if remaining_weight <= 0:
                step.append((box, lifter_indices))
                boxes.remove(box)
        
        if step:
            steps.append(step)
        if not boxes:
            break
    
    return steps if not boxes else None

boxes = [55, 244, 173, 293, 90, 126, 340, 250, 66, 143, 103, 244, 76, 166, 130, 216, 54, 196, 245, 307]
lifters = [142, 178, 196, 52, 101, 144, 50]

steps = assign_lifters_to_boxes(boxes, lifters)

if steps:
    output = ""
    for i, step in enumerate(steps):
        output += f"Step {i + 1}: {step}\n"
    print(f"<<<{output}>>>")
else:
    print("<<<No valid solution found within 6 steps>>>")