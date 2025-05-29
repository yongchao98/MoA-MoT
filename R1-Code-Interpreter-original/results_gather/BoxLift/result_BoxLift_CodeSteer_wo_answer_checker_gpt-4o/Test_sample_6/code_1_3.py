from functools import lru_cache

boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifters = [61, 45, 66, 38, 52]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def find_steps(remaining_boxes, used_lifters):
    if not remaining_boxes:
        return []

    for i, box in enumerate(remaining_boxes):
        total_capacity = 0
        lifter_indices = []
        
        for j, lifter in enumerate(lifters):
            if not used_lifters[j] and total_capacity + lifter <= box:
                total_capacity += lifter
                lifter_indices.append(j)
                if total_capacity == box:
                    break
        
        if total_capacity >= box:
            new_used_lifters = list(used_lifters)
            for idx in lifter_indices:
                new_used_lifters[idx] = True
            result = find_steps(remaining_boxes[:i] + remaining_boxes[i+1:], tuple(new_used_lifters))
            if result is not None:
                return [(box, lifter_indices)] + result

    return None

steps = []
remaining_boxes = tuple(boxes)
used_lifters = tuple([False] * len(lifters))

while remaining_boxes:
    step = find_steps(remaining_boxes, used_lifters)
    if step is None:
        break
    steps.append(step)
    remaining_boxes = tuple(box for box, _ in step if box not in remaining_boxes)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")