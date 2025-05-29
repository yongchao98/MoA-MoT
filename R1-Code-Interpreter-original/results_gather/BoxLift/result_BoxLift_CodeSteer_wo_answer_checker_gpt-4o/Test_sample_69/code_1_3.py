from functools import lru_cache

boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

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
            new_used_lifters = used_lifters
            for idx in lifter_indices:
                new_used_lifters |= (1 << idx)
            remaining_boxes_copy = remaining_boxes[:i] + remaining_boxes[i+1:]
            result = find_steps(tuple(remaining_boxes_copy), new_used_lifters)
            if result is not None:
                return [(box, lifter_indices)] + result

    return None

steps = []
while boxes:
    step = find_steps(tuple(boxes), 0)
    if step is None:
        break
    steps.append(step)
    for box, _ in step:
        boxes.remove(box)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")