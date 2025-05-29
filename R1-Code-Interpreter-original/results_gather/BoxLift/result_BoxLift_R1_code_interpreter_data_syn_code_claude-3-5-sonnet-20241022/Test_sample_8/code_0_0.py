from itertools import combinations
import copy

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                return list(combo)
    return None

# Initialize data
boxes = [362, 202, 192, 149, 357, 335, 209, 193, 334, 375, 374, 167, 276, 70, 282, 201, 242, 158, 236, 212]
lifters = [174, 130, 182, 186, 80, 149, 87]

# Sort boxes in descending order
boxes = sorted(boxes, reverse=True)
result = []
remaining_boxes = boxes.copy()

while remaining_boxes:
    step = []
    available_lifters = list(range(len(lifters)))
    current_boxes = remaining_boxes.copy()
    
    # Try to assign boxes in current step
    for box in current_boxes:
        if not available_lifters:
            break
            
        # Get lifter combination that can lift this box
        lifter_combo = can_lift([lifters[i] for i in available_lifters], box)
        
        if lifter_combo:
            # Convert relative indices to absolute indices
            actual_lifters = [available_lifters[i] for i in lifter_combo]
            step.append((box, actual_lifters))
            # Remove used lifters
            for idx in sorted(lifter_combo, reverse=True):
                available_lifters.pop(idx)
            remaining_boxes.remove(box)
    
    if step:
        result.append(step)

# Format output
output = ""
for i, step in enumerate(result, 1):
    output += f"Step {i}: {step}\n"

print(output.strip())