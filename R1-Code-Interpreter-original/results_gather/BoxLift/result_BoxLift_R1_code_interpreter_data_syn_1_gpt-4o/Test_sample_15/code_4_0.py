from itertools import combinations

boxes = [50, 70, 83, 17, 21, 27, 53, 31]
lifters = [42, 53, 74, 78]

# Sort boxes in descending order
boxes.sort(reverse=True)

def can_lift(box_weight, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box_weight

def find_lifter_combinations(box_weight, available_lifters):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(available_lifters, r):
            if can_lift(box_weight, combo):
                return combo
    return None

def assign_boxes_to_steps(boxes, lifters, max_steps):
    steps = []
    remaining_boxes = boxes[:]
    while remaining_boxes and len(steps) < max_steps:
        step = []
        available_lifters = list(range(len(lifters)))
        for box in remaining_boxes[:]:
            lifter_combo = find_lifter_combinations(box, available_lifters)
            if lifter_combo:
                step.append((box, list(lifter_combo)))
                for index in lifter_combo:
                    available_lifters.remove(index)
                remaining_boxes.remove(box)
        steps.append(step)
    return steps if not remaining_boxes else "Not possible to lift all boxes in 3 steps or less"

steps = assign_boxes_to_steps(boxes, lifters, 3)
print(steps)