from itertools import combinations

def can_lift(lifters, weight):
    # Check all possible combinations of lifters
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= weight:
                return list(combo)
    return None

def solve_box_lifting(boxes, lifter_capacities):
    # Sort boxes from heaviest to lightest
    boxes = sorted(boxes, reverse=True)
    steps = []
    
    while boxes:
        step = []
        available_lifters = lifter_capacities.copy()
        used_lifter_indices = set()
        
        # Try to assign as many boxes as possible in current step
        i = 0
        while i < len(boxes):
            box = boxes[i]
            
            # Find available lifters that can lift this box
            lifter_combo = None
            best_combo = None
            min_waste = float('inf')
            
            # Try all possible combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                for combo_indices in combinations(range(len(available_lifters)), r):
                    # Skip if any lifter in combo is already used
                    if any(idx in used_lifter_indices for idx in combo_indices):
                        continue
                    
                    combo_capacity = sum(available_lifters[j] for j in combo_indices)
                    if combo_capacity >= box and combo_capacity - box < min_waste:
                        min_waste = combo_capacity - box
                        best_combo = combo_indices
                        
                        # If waste is 0 or very small, use this combination
                        if min_waste <= 10:
                            break
                
                if best_combo and min_waste <= 10:
                    break
            
            if best_combo:
                # Add assignment to step
                step.append((box, [idx for idx in best_combo]))
                # Mark lifters as used
                used_lifter_indices.update(best_combo)
                # Remove box from list
                boxes.pop(i)
            else:
                i += 1
        
        if step:
            steps.append(step)
    
    # Format output
    result = ""
    for i, step in enumerate(steps, 1):
        result += f"Step {i}: {step}\n"
    
    return result.strip()

# Input data
boxes = [151, 393, 143, 179, 88, 173, 344, 212, 308, 258, 113, 261, 371, 335, 236, 258, 124, 262, 108, 209, 232, 128, 266, 232]
lifters = [84, 94, 87, 108, 198, 95, 121]

print(solve_box_lifting(boxes, lifters))