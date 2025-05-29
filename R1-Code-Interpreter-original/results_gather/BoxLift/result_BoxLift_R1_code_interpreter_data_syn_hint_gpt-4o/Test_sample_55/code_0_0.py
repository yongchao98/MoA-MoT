boxes = [24, 78, 95, 23, 57, 36, 40, 26]
lifters = [56, 51, 66, 65]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

steps = []
used_boxes = [False] * len(boxes)

def can_lift(box_weight, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box_weight

def find_combination(box_weight, used_lifters):
    n = len(lifters)
    for i in range(1, 1 << n):
        combination = [j for j in range(n) if (i & (1 << j)) and not used_lifters[j]]
        if can_lift(box_weight, combination):
            return combination
    return None

for _ in range(3):  # Maximum 3 steps
    step = []
    used_lifters = [False] * len(lifters)
    for i, box in enumerate(boxes):
        if not used_boxes[i]:
            combination = find_combination(box, used_lifters)
            if combination is not None:
                step.append((box, combination))
                used_boxes[i] = True
                for idx in combination:
                    used_lifters[idx] = True
    if step:
        steps.append(step)
    if all(used_boxes):
        break

print(steps)