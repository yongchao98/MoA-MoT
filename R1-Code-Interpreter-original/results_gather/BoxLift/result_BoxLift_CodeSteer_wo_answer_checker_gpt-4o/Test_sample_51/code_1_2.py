from functools import lru_cache

boxes = [142, 109, 62, 158, 37, 50, 180, 128, 61, 123, 44, 51]
lifters = [93, 71, 82, 91]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def can_lift(box_index, used_lifters):
    if box_index == len(boxes):
        return True, []

    box = boxes[box_index]
    lifter_indices = []
    total_capacity = 0

    for i, capacity in enumerate(lifters):
        if not (used_lifters & (1 << i)) and total_capacity + capacity <= box:
            total_capacity += capacity
            lifter_indices.append(i)
            if total_capacity >= box:
                break

    if total_capacity >= box:
        new_used_lifters = used_lifters
        for idx in lifter_indices:
            new_used_lifters |= (1 << idx)
        can_lift_next, next_steps = can_lift(box_index + 1, new_used_lifters)
        if can_lift_next:
            return True, [(box, lifter_indices)] + next_steps

    return False, []

steps = []
used_lifters = 0
for step in range(5):
    can_lift_all, step_result = can_lift(0, used_lifters)
    if can_lift_all:
        steps.append(step_result)
        break
    else:
        steps.append(step_result)
        for _, lifter_indices in step_result:
            for idx in lifter_indices:
                used_lifters |= (1 << idx)

# Print the steps
output = ""
for i, step in enumerate(steps):
    output += f"Step {i + 1}: {step}\n"

print(f"<<<{output.strip()}>>>")