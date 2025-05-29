from itertools import combinations

def find_solution():
    boxes = [100, 22, 67, 79, 32, 31, 67, 37]
    lifter_caps = [58, 78, 68, 75]
    box_indices = list(range(len(boxes)))
    n_lifters = len(lifter_caps)
    
    def try_assign_step(remaining_boxes, used_steps):
        if not remaining_boxes:  # All boxes assigned
            return used_steps
        if len(used_steps) >= 3:  # Exceeded max steps
            return None
            
        # Sort remaining boxes by weight (descending)
        boxes_list = sorted(list(remaining_boxes), key=lambda x: boxes[x], reverse=True)
        current_step = []
        used_lifters = set()
        new_remaining = remaining_boxes.copy()
        
        # Try to assign as many boxes as possible in this step
        for box_idx in boxes_list:
            if len(used_lifters) >= n_lifters:
                break
                
            box_weight = boxes[box_idx]
            available_lifters = [i for i in range(n_lifters) if i not in used_lifters]
            
            # Try different combinations of lifters
            for r in range(1, len(available_lifters) + 1):
                for lifter_combo in combinations(available_lifters, r):
                    combo_capacity = sum(lifter_caps[i] for i in lifter_combo)
                    if combo_capacity >= box_weight:
                        current_step.append((box_weight, list(lifter_combo)))
                        used_lifters.update(lifter_combo)
                        new_remaining.remove(box_idx)
                        break
                else:
                    continue
                break
        
        if current_step:
            next_solution = try_assign_step(new_remaining, used_steps + [current_step])
            if next_solution is not None:
                return next_solution
        
        return None

    # Start with all boxes
    initial_boxes = set(box_indices)
    solution = try_assign_step(initial_boxes, [])
    
    if solution:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
        return True
    return False

find_solution()