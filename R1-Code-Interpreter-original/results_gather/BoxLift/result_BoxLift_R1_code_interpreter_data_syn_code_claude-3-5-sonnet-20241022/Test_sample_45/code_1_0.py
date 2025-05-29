from itertools import combinations

def can_lift(lifters_cap, weight):
    return sum(lifters_cap) >= weight

def find_solution():
    boxes = [100, 22, 67, 79, 32, 31, 67, 37]
    lifter_caps = [58, 78, 68, 75]
    n_lifters = len(lifter_caps)
    
    # Track boxes by their indices
    box_indices = list(range(len(boxes)))
    # Sort boxes by weight (descending)
    box_indices.sort(key=lambda x: boxes[x], reverse=True)
    
    solution = []
    remaining_boxes = set(box_indices)
    
    while remaining_boxes and len(solution) < 3:
        step = []
        used_lifters = set()
        
        # First handle heavy boxes that need multiple lifters
        for box_idx in list(remaining_boxes):
            if len(used_lifters) >= n_lifters:
                break
                
            box_weight = boxes[box_idx]
            available_lifters = sorted(
                [i for i in range(n_lifters) if i not in used_lifters],
                key=lambda x: lifter_caps[x],
                reverse=True
            )
            
            # Try combinations of lifters, starting with minimum needed
            for r in range(1, len(available_lifters) + 1):
                found = False
                for lifter_combo in combinations(available_lifters, r):
                    combo_capacity = sum(lifter_caps[i] for i in lifter_combo)
                    if combo_capacity >= box_weight:
                        step.append((box_weight, list(lifter_combo)))
                        used_lifters.update(lifter_combo)
                        remaining_boxes.remove(box_idx)
                        found = True
                        break
                if found:
                    break
        
        if step:
            solution.append(step)
        else:
            # If we can't make progress, break
            return False
    
    if not remaining_boxes:
        # Format and print the solution
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

find_solution()