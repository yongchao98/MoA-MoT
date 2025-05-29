from itertools import combinations
import copy

def can_lift(lifters, weight):
    if not lifters:
        return False
    return sum(lifters) >= weight

def find_lifter_combination(available_lifters, weight):
    for r in range(1, len(available_lifters) + 1):
        for combo in combinations(range(len(available_lifters)), r):
            lifter_weights = [available_lifters[i] for i in combo]
            if sum(lifter_weights) >= weight:
                return list(combo)
    return None

boxes = [110, 230, 255, 47, 133, 280, 271, 275, 155, 80, 169, 89, 299, 241, 187, 234]
lifters = [160, 110, 114, 55, 54, 60]

# Sort boxes in descending order
boxes_with_index = list(enumerate(boxes))
boxes_with_index.sort(key=lambda x: x[1], reverse=True)

steps = []
remaining_boxes = boxes_with_index.copy()

while remaining_boxes:
    step = []
    available_lifters = list(range(len(lifters)))
    
    current_boxes = remaining_boxes.copy()
    for box_idx, box_weight in current_boxes:
        if not available_lifters:
            break
            
        lifter_combo = find_lifter_combination([lifters[i] for i in available_lifters], box_weight)
        if lifter_combo:
            # Convert relative indices to absolute indices
            absolute_lifter_indices = [available_lifters[i] for i in lifter_combo]
            step.append((box_weight, absolute_lifter_indices))
            # Remove used lifters
            for idx in sorted(lifter_combo, reverse=True):
                available_lifters.pop(idx)
            remaining_boxes.remove((box_idx, box_weight))
    
    steps.append(step)

# Format output
result = ""
for i, step in enumerate(steps, 1):
    result += f"Step {i}: {step}\n"

print(result.strip())