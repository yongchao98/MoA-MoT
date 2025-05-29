from functools import lru_cache

boxes = [281, 159, 239, 61, 219, 44, 225, 192, 74, 369, 318, 366, 62, 186, 373, 243, 131, 172, 66, 167, 400, 330, 351, 245]
lifters = [164, 182, 169, 127, 118, 60]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def find_steps(remaining_boxes, used_lifters):
    if not remaining_boxes:
        return []

    for i, box in enumerate(remaining_boxes):
        lifter_indices = []
        total_capacity = 0
        for j, lifter in enumerate(lifters):
            if not (used_lifters & (1 << j)) and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity >= box:
                    break

        if total_capacity >= box:
            next_used_lifters = used_lifters
            for idx in lifter_indices:
                next_used_lifters |= (1 << idx)
            next_steps = find_steps(remaining_boxes[:i] + remaining_boxes[i+1:], next_used_lifters)
            if next_steps is not None:
                return [(box, lifter_indices)] + next_steps

    return None

steps = []
for step in range(9):
    step_lifts = find_steps(tuple(boxes), 0)
    if step_lifts is None:
        break
    steps.append(step_lifts)
    for box, lifter_indices in step_lifts:
        for idx in lifter_indices:
            lifters[idx] = 0  # Mark these lifters as used for this step
    boxes = [box for box, _ in step_lifts if box not in boxes]

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")