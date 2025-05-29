from itertools import combinations

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_lifter_combinations(lifter_capacities, weight, used_lifters, max_lifters=5):
    available_lifters = [i for i in range(len(lifter_capacities)) if i not in used_lifters]
    for r in range(1, min(len(available_lifters) + 1, max_lifters + 1)):
        for comb in combinations(available_lifters, r):
            if sum(lifter_capacities[i] for i in comb) >= weight:
                return list(comb)
    return None

def solve_box_lifting(boxes, lifter_capacities):
    # Sort boxes in descending order
    boxes = sorted(boxes, reverse=True)
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes and len(steps) < 10:
        step = []
        used_lifters = set()
        
        # Try to assign boxes to lifters
        i = 0
        while i < len(remaining_boxes):
            box = remaining_boxes[i]
            lifters = find_lifter_combinations(lifter_capacities, box, used_lifters)
            
            if lifters:
                step.append((box, lifters))
                used_lifters.update(lifters)
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format output
    result = []
    for i, step in enumerate(steps, 1):
        result.append(f"Step {i}: {step}")
    
    return "\n".join(result)

# Input data
boxes = [44, 178, 94, 176, 159, 199, 197, 133, 173, 186, 196, 80]
lifter_capacities = [61, 45, 66, 38, 52]

print(solve_box_lifting(boxes, lifter_capacities))