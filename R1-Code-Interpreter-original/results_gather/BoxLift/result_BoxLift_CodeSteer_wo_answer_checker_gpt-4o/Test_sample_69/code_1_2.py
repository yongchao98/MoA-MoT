from functools import lru_cache

boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

@lru_cache(None)
def can_lift(box_index, used_lifters):
    if box_index >= len(boxes):
        return True, []

    box = boxes[box_index]
    for i in range(len(lifters)):
        if not (used_lifters & (1 << i)) and lifters[i] >= box:
            result, steps = can_lift(box_index + 1, used_lifters | (1 << i))
            if result:
                return True, [(box, [i])] + steps

    for i in range(len(lifters)):
        if not (used_lifters & (1 << i)):
            for j in range(i + 1, len(lifters)):
                if not (used_lifters & (1 << j)) and lifters[i] + lifters[j] >= box:
                    result, steps = can_lift(box_index + 1, used_lifters | (1 << i) | (1 << j))
                    if result:
                        return True, [(box, [i, j])] + steps

    return False, []

steps = []
current_step = 0
while boxes:
    result, step = can_lift(0, 0)
    if not result:
        break
    steps.append(step)
    for box, lifter_indices in step:
        boxes.remove(box)
    current_step += 1

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i + 1}: {step}")

print("<<<" + "\n".join(output) + ">>>")