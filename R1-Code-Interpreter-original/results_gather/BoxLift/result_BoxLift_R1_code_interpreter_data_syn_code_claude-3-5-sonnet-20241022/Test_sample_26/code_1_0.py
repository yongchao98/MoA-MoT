def solve_box_lifting():
    boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
    lifters = [57, 61, 104, 98, 70]
    
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    
    # Pre-calculate two strongest lifters combined capacity
    strongest_pair = lifters[2] + lifters[3]  # 104 + 98 = 202
    
    solution = []
    remaining_boxes = boxes[:]
    
    while remaining_boxes:
        step = []
        used_lifters = set()
        
        # Process remaining boxes in current step
        i = 0
        while i < len(remaining_boxes):
            box = remaining_boxes[i]
            
            # For very heavy boxes (>140), use strongest pair
            if box > 140 and 2 not in used_lifters and 3 not in used_lifters:
                if strongest_pair >= box:
                    step.append((box, [2, 3]))
                    used_lifters.add(2)
                    used_lifters.add(3)
                    remaining_boxes.pop(i)
                    continue
            
            # For medium boxes, try single lifters
            for j in range(len(lifters)):
                if j not in used_lifters and lifters[j] >= box:
                    step.append((box, [j]))
                    used_lifters.add(j)
                    remaining_boxes.pop(i)
                    break
            else:
                # If no single lifter can handle it, try pairs
                for j in range(len(lifters)):
                    if j in used_lifters:
                        continue
                    for k in range(j + 1, len(lifters)):
                        if k in used_lifters:
                            continue
                        if lifters[j] + lifters[k] >= box:
                            step.append((box, [j, k]))
                            used_lifters.add(j)
                            used_lifters.add(k)
                            remaining_boxes.pop(i)
                            break
                    else:
                        continue
                    break
            i += 1
            
        if step:
            solution.append(step)
        
        if len(solution) > 4:
            return "No solution in 4 steps"
    
    # Format solution
    result = ""
    for i, step in enumerate(solution, 1):
        result += f"Step {i}: {step}\n"
    return result.strip()

print(solve_box_lifting())