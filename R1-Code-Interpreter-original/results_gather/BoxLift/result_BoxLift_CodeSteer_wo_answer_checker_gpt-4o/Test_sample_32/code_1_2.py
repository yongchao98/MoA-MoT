def can_lift(box, lifters, used_lifters):
    # Try to find a combination of lifters that can lift the box
    def backtrack(remaining_weight, start, current_lifters):
        if remaining_weight <= 0:
            return current_lifters
        for i in range(start, len(lifters)):
            if not used_lifters[i] and lifters[i] <= remaining_weight:
                used_lifters[i] = True
                result = backtrack(remaining_weight - lifters[i], i + 1, current_lifters + [i])
                if result is not None:
                    return result
                used_lifters[i] = False
        return None

    return backtrack(box, 0, [])

def assign_lifters_to_boxes(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = []

    while boxes:
        step = []
        used_lifters = [False] * len(lifters)
        
        for box in boxes[:]:
            lifter_indices = can_lift(box, lifters, used_lifters)
            if lifter_indices is not None:
                step.append((box, lifter_indices))
                for index in lifter_indices:
                    used_lifters[index] = True
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