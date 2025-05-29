def solve_box_lifting():
    boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
    lifters = [81, 98, 62, 46]
    
    # Sort boxes by weight (descending)
    box_indices = list(range(len(boxes)))
    box_indices.sort(key=lambda x: boxes[x], reverse=True)
    
    # Initialize solution
    solution = []
    unassigned_boxes = set(box_indices)
    
    while unassigned_boxes and len(solution) < 7:
        step_assignments = []
        available_lifters = set(range(len(lifters)))
        
        # Process remaining boxes
        for box_idx in list(unassigned_boxes):
            box_weight = boxes[box_idx]
            
            # Skip if no lifters available
            if not available_lifters:
                continue
                
            # Find best lifter combination
            best_lifters = None
            min_waste = float('inf')
            
            # Try all possible combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                if r > 3:  # Limit combinations to maximum 3 lifters per box
                    break
                    
                lifter_list = sorted(list(available_lifters), 
                                   key=lambda x: lifters[x], reverse=True)
                
                for i in range(len(lifter_list) - r + 1):
                    combo = lifter_list[i:i+r]
                    total_capacity = sum(lifters[j] for j in combo)
                    
                    if total_capacity >= box_weight:
                        waste = total_capacity - box_weight
                        if waste < min_waste:
                            min_waste = waste
                            best_lifters = combo
                            
                        if waste == 0:  # Perfect match
                            break
                            
                if best_lifters and min_waste == 0:
                    break
            
            if best_lifters:
                step_assignments.append((box_weight, list(best_lifters)))
                available_lifters -= set(best_lifters)
                unassigned_boxes.remove(box_idx)
        
        if step_assignments:
            solution.append(f"Step {len(solution) + 1}: {step_assignments}")
        else:
            break
    
    if not unassigned_boxes:
        print("<<<" + "\n".join(solution) + ">>>")
    else:
        print("No solution found within 7 steps")

solve_box_lifting()