def find_lifter_combinations(weight: int, lifters: list, available_lifters: list) -> list:
    valid_combinations = []
    available_indices = [i for i, is_available in enumerate(available_lifters) if is_available]
    
    for r in range(1, len(available_indices) + 1):
        for combo in itertools.combinations(available_indices, r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combinations.append(list(combo))
    return sorted(valid_combinations, key=len)  # Sort by number of lifters needed

def solve_box_lifting(boxes: list, lifters: list) -> list:
    # Sort boxes in descending order
    boxes = sorted([(weight, i) for i, weight in enumerate(boxes)], reverse=True)
    solution = []
    unassigned_boxes = boxes.copy()
    
    while unassigned_boxes:
        step = []
        available_lifters = [True] * len(lifters)
        
        # Try to assign as many boxes as possible in current step
        remaining_boxes = unassigned_boxes.copy()
        for box_weight, box_idx in remaining_boxes:
            # Skip if no lifters available
            if not any(available_lifters):
                break
                
            # Find valid lifter combinations for this box
            valid_combinations = find_lifter_combinations(box_weight, lifters, available_lifters)
            
            if valid_combinations:
                # Use the combination with fewest lifters
                best_combo = valid_combinations[0]
                
                # Add assignment to current step
                step.append((box_weight, best_combo))
                
                # Mark lifters as unavailable
                for lifter_idx in best_combo:
                    available_lifters[lifter_idx] = False
                    
                # Remove box from unassigned
                unassigned_boxes.remove((box_weight, box_idx))
        
        if step:
            solution.append(step)
        
        # Safety check to prevent infinite loop
        if len(solution) > 7:
            return None
    
    return solution

# Input data
boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Solve and format output
solution = solve_box_lifting(boxes, lifters)

if solution is None:
    print("<<<No solution found within 7 steps>>>")
else:
    output = []
    for step_idx, step in enumerate(solution, 1):
        output.append(f"Step {step_idx}: {step}")
    print("<<<" + "\n".join(output) + ">>>")