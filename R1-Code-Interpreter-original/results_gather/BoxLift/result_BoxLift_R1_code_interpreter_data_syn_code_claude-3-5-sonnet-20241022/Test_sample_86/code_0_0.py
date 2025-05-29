from itertools import combinations
import copy

def can_lift(box_weight, lifters_combination, capacities):
    return sum(capacities[i] for i in lifters_combination) >= box_weight

def find_solution():
    boxes = [57, 112, 121, 49, 166, 143, 169, 43, 90, 66, 23, 143]
    capacities = [103, 49, 95, 97, 98]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    boxes = [(i, w) for i, w in boxes]  # Keep original index
    
    steps = []
    while boxes:
        step = []
        available_lifters = set(range(len(capacities)))
        
        # Try to assign as many boxes as possible in this step
        remaining_boxes = []
        for box_idx, box_weight in boxes:
            if not available_lifters:
                remaining_boxes.append((box_idx, box_weight))
                continue
                
            # Try different combinations of available lifters
            assigned = False
            for size in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifters in combinations(available_lifters, size):
                    if can_lift(box_weight, lifters, capacities):
                        step.append((box_weight, list(lifters)))
                        available_lifters -= set(lifters)
                        assigned = True
                        break
            
            if not assigned:
                remaining_boxes.append((box_idx, box_weight))
        
        steps.append(step)
        boxes = remaining_boxes
        
        if len(steps) > 4:  # Check if we exceed maximum steps
            return None
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

solution = find_solution()
print(solution)