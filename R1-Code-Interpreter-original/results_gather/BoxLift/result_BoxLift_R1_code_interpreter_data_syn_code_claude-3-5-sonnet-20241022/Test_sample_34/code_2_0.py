import itertools
from collections import defaultdict

def find_lifting_solution():
    boxes = [399, 277, 278, 282, 209, 309, 281, 53, 325, 239, 363, 363, 308, 244, 392, 323, 255, 79, 280, 336, 395, 318, 203, 364]
    lifters = [125, 145, 57, 73, 81, 70]
    boxes.sort(reverse=True)
    
    # Precompute all valid lifter combinations and their total strength
    lifter_combos = {}
    for r in range(1, len(lifters) + 1):
        for combo in itertools.combinations(range(len(lifters)), r):
            lifter_combos[combo] = sum(lifters[i] for i in combo)
    
    def find_optimal_assignment(available_boxes, available_lifters):
        assignments = []
        used_lifters = set()
        
        # Try to assign heavy boxes first
        for box in available_boxes[:]:
            if not available_lifters - used_lifters:
                break
                
            # Find the most efficient lifter combination for this box
            best_combo = None
            min_excess = float('inf')
            for combo in itertools.combinations(available_lifters - used_lifters, min(4, len(available_lifters - used_lifters))):
                if sum(lifters[i] for i in combo) >= box:
                    excess = sum(lifters[i] for i in combo) - box
                    if excess < min_excess:
                        min_excess = excess
                        best_combo = combo
            
            if best_combo:
                assignments.append((box, list(best_combo)))
                used_lifters.update(best_combo)
                available_boxes.remove(box)
        
        return assignments if assignments else None
    
    steps = []
    remaining_boxes = boxes.copy()
    
    while remaining_boxes:
        available_lifters = set(range(len(lifters)))
        step_assignments = find_optimal_assignment(remaining_boxes, available_lifters)
        
        if step_assignments:
            steps.append(step_assignments)
        else:
            # If no assignment found, something went wrong
            break
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

find_lifting_solution()