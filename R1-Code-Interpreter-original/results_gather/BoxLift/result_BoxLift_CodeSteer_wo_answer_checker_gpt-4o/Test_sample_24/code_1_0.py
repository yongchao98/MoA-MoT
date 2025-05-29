from itertools import combinations

boxes = [119, 110, 359, 148, 237, 71, 309, 213, 117, 394, 79, 353, 230, 77, 242, 239, 66, 259, 212, 358, 194, 397, 274, 101]
lifters = [60, 103, 156, 108, 179, 110]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Precompute possible combinations of lifters for different weight ranges
def precompute_lifter_combinations(lifters):
    lifter_combinations = {}
    for r in range(1, len(lifters) + 1):
        for combo in combinations(enumerate(lifters), r):
            indices, capacities = zip(*combo)
            total_capacity = sum(capacities)
            if total_capacity not in lifter_combinations:
                lifter_combinations[total_capacity] = []
            lifter_combinations[total_capacity].append(list(indices))
    return lifter_combinations

lifter_combinations = precompute_lifter_combinations(lifters)

def can_lift(box_weight, lifters, used):
    for capacity, combos in sorted(lifter_combinations.items(), reverse=True):
        if capacity >= box_weight:
            for combo in combos:
                if not any(used[i] for i in combo):
                    return combo
    return None

def branch_and_bound(boxes, lifters, step, steps, used, memo):
    if step >= 9:
        return False
    if not boxes:
        return True
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    current_step = []
    remaining_boxes = boxes[:]
    
    for box in boxes:
        lifter_indices = can_lift(box, lifters, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            for i in lifter_indices:
                used[i] = True
            remaining_boxes.remove(box)
    
    if current_step:
        steps.append(current_step)
        if branch_and_bound(remaining_boxes, lifters, step + 1, steps, [False] * len(lifters), memo):
            memo[state] = True
            return True
        steps.pop()
    
    memo[state] = False
    return False

steps = []
used = [False] * len(lifters)
memo = {}
if branch_and_bound(boxes, lifters, 0, steps, used, memo):
    output = "\n".join(f"Step {i+1}: {step}" for i, step in enumerate(steps))
    print(f"<<<{output}>>>")
else:
    print("No solution within 9 steps.")