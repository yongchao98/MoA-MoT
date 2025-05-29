boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
lifters = [167, 138, 133, 87, 142, 167, 140]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box_index, used_lifters, steps, current_step):
    if box_index == len(boxes):
        return steps <= 7, current_step
    if steps > 7:
        return False, []

    box_weight = boxes[box_index]
    lifter_indices = []
    total_capacity = 0

    for i in range(len(lifters)):
        if not used_lifters[i] and total_capacity < box_weight:
            lifter_indices.append(i)
            total_capacity += lifters[i]
            used_lifters[i] = True
        if total_capacity >= box_weight:
            break

    if total_capacity >= box_weight:
        current_step.append((box_weight, lifter_indices))
        for i in lifter_indices:
            used_lifters[i] = False
        success, result = can_lift(box_index + 1, used_lifters, steps + 1, [])
        if success:
            return True, [current_step] + result

    return False, []

def find_solution():
    used_lifters = [False] * len(lifters)
    success, steps = can_lift(0, used_lifters, 0, [])
    if success:
        output = "<<<"
        for i, step in enumerate(steps):
            output += f"Step {i + 1}: {step}\n"
        output += ">>>"
        return output
    else:
        return "It's not possible to lift all boxes in 7 or fewer steps."

print(find_solution())