from functools import lru_cache

boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def can_lift(box_index, used_lifters):
    if box_index == len(boxes):
        return True, []

    box_weight = boxes[box_index]
    lifter_indices = []
    total_capacity = 0

    for i, lifter in enumerate(lifters):
        if not (used_lifters & (1 << i)) and total_capacity + lifter <= box_weight:
            total_capacity += lifter
            lifter_indices.append(i)
            if total_capacity >= box_weight:
                break

    if total_capacity >= box_weight:
        next_used_lifters = used_lifters
        for idx in lifter_indices:
            next_used_lifters |= (1 << idx)
        can_lift_next, next_steps = can_lift(box_index + 1, next_used_lifters)
        if can_lift_next:
            return True, [(box_weight, lifter_indices)] + next_steps

    return False, []

steps = []
used_lifters = 0
for step in range(9):
    can_lift_all, step_lifts = can_lift(0, used_lifters)
    if can_lift_all:
        steps.append(step_lifts)
        break
    else:
        steps.append(step_lifts)
        for _, lifter_indices in step_lifts:
            for idx in lifter_indices:
                used_lifters |= (1 << idx)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")