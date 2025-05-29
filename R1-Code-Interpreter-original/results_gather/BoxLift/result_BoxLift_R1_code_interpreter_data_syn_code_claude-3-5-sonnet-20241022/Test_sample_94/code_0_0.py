from itertools import combinations
import copy

# Input data
boxes = [187, 344, 112, 66, 384, 247, 184, 298, 131, 51, 62, 255, 120, 357, 399, 287, 231, 161, 336, 256, 328, 239, 365, 245]
lifters = [165, 85, 52, 82, 114, 75]

def get_lifter_combinations(lifters):
    combinations_dict = {}
    # Get all possible combinations of lifters and their total strength
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            total_strength = sum(lifters[i] for i in combo)
            combinations_dict[combo] = total_strength
    return combinations_dict

def solve_box_lifting():
    # Sort boxes in descending order
    remaining_boxes = sorted(boxes, reverse=True)
    steps = []
    lifter_combinations = get_lifter_combinations(lifters)
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes to available lifters
        for box in remaining_boxes[:]:
            if not available_lifters:
                break
                
            # Find the best combination of available lifters for this box
            best_combo = None
            min_lifters = float('inf')
            
            for combo, strength in lifter_combinations.items():
                if (strength >= box and 
                    len(combo) < min_lifters and 
                    all(l in available_lifters for l in combo)):
                    best_combo = combo
                    min_lifters = len(combo)
            
            if best_combo:
                step.append((box, list(best_combo)))
                for lifter in best_combo:
                    available_lifters.remove(lifter)
                remaining_boxes.remove(box)
        
        if step:
            steps.append(step)
    
    # Format the output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

print(solve_box_lifting())