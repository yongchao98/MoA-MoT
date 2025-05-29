import itertools

def solve_box_lifting():
    boxes = [295, 145, 158, 66, 383, 291, 156, 366, 271, 121, 289, 388, 276, 247, 68, 248, 167, 252, 126, 159]
    lifters = [144, 152, 155, 109, 115, 152, 162]
    
    # Sort boxes in descending order
    boxes = sorted(enumerate(boxes), key=lambda x: x[1], reverse=True)
    remaining_boxes = boxes.copy()
    steps = []
    
    def get_lifter_combinations():
        combinations = []
        for r in range(1, len(lifters) + 1):
            combinations.extend(list(itertools.combinations(range(len(lifters)), r)))
        return combinations
    
    def can_lift(weight, lifter_combo):
        return sum(lifters[i] for i in lifter_combo) >= weight
    
    while remaining_boxes and len(steps) < 6:
        step = []
        used_lifters = set()
        
        # Try to assign boxes in current step
        boxes_to_remove = []
        for box_idx, box_weight in remaining_boxes:
            if not (set(range(len(lifters))) - used_lifters):
                break
                
            # Get all possible lifter combinations that haven't been used in this step
            available_lifters = list(set(range(len(lifters))) - used_lifters)
            lifter_combos = []
            for r in range(1, len(available_lifters) + 1):
                lifter_combos.extend(list(itertools.combinations(available_lifters, r)))
            
            # Find the smallest combination that can lift the box
            valid_combos = [combo for combo in lifter_combos if can_lift(box_weight, combo)]
            if valid_combos:
                best_combo = min(valid_combos, key=len)
                step.append((box_weight, list(best_combo)))
                used_lifters.update(best_combo)
                boxes_to_remove.append((box_idx, box_weight))
        
        # Remove assigned boxes
        for box in boxes_to_remove:
            remaining_boxes.remove(box)
        
        steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    print(result.strip())

solve_box_lifting()