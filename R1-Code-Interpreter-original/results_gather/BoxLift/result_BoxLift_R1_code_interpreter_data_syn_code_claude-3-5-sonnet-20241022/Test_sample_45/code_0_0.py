from itertools import combinations

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_solution():
    boxes = [100, 22, 67, 79, 32, 31, 67, 37]
    lifter_caps = [58, 78, 68, 75]
    n_lifters = len(lifter_caps)
    
    # Sort boxes in descending order for better assignment
    box_indices = list(range(len(boxes)))
    box_indices.sort(key=lambda x: boxes[x], reverse=True)
    
    solution = []
    remaining_boxes = set(box_indices)
    
    while remaining_boxes and len(solution) < 3:
        step = []
        used_lifters = set()
        
        # Try to assign boxes in current step
        for box_idx in list(remaining_boxes):
            if len(used_lifters) == n_lifters:
                break
                
            box_weight = boxes[box_idx]
            # Try different combinations of available lifters
            available_lifters = [i for i in range(n_lifters) if i not in used_lifters]
            
            for r in range(1, len(available_lifters) + 1):
                for lifter_combo in combinations(available_lifters, r):
                    lifter_sum = sum(lifter_caps[i] for i in lifter_combo)
                    if lifter_sum >= box_weight:
                        step.append((box_weight, list(lifter_combo)))
                        used_lifters.update(lifter_combo)
                        remaining_boxes.remove(box_idx)
                        break
                else:
                    continue
                break
        
        if step:
            solution.append(step)
        else:
            break
    
    if not remaining_boxes:
        # Format the solution
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

find_solution()