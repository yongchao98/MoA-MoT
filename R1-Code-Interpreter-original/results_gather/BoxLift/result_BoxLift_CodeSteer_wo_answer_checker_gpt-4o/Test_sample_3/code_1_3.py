from functools import lru_cache

boxes = [159, 197, 211, 302, 51, 251, 299, 229, 384, 196, 377, 283, 127, 312, 228, 280, 153, 199, 169, 151]
lifters = [52, 80, 175, 75, 136, 143]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def can_lift(remaining_boxes, used_lifters_mask):
    if not remaining_boxes:
        return []

    for i, box in enumerate(remaining_boxes):
        total_capacity = 0
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not (used_lifters_mask & (1 << j)) and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity >= box:
                    break
        
        if total_capacity >= box:
            new_used_lifters_mask = used_lifters_mask
            for idx in lifter_indices:
                new_used_lifters_mask |= (1 << idx)
            result = can_lift(remaining_boxes[:i] + remaining_boxes[i+1:], new_used_lifters_mask)
            if result is not None:
                return [(box, lifter_indices)] + result

    return None

steps = []
used_lifters_mask = 0
while boxes:
    step = can_lift(tuple(boxes), used_lifters_mask)
    if step is None:
        break
    steps.append(step)
    for box, lifter_indices in step:
        boxes.remove(box)
        for idx in lifter_indices:
            used_lifters_mask |= (1 << idx)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")