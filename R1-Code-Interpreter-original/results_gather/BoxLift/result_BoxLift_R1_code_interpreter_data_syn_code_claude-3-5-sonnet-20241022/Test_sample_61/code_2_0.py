from itertools import combinations

def get_all_lifter_combinations(lifters):
    all_combos = {}
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            strength = sum(lifters[i] for i in combo)
            if strength not in all_combos or len(combo) < len(all_combos[strength]):
                all_combos[strength] = list(combo)
    return all_combos

def find_best_combination(weight, lifter_combos, available_lifters):
    best_combo = None
    min_lifters = float('inf')
    
    for strength, combo in lifter_combos.items():
        if strength >= weight and all(l in available_lifters for l in combo):
            if len(combo) < min_lifters:
                min_lifters = len(combo)
                best_combo = combo
    
    return best_combo

def solve_box_lifting():
    boxes = [128, 144, 42, 49, 132, 118, 32, 48, 21, 94, 85, 51]
    lifters = [85, 90, 40, 102]
    
    # Pre-calculate all valid lifter combinations and their total strengths
    lifter_combinations = get_all_lifter_combinations(lifters)
    
    # Sort boxes by weight (heaviest first)
    boxes = sorted(enumerate(boxes), key=lambda x: (-x[1], x[0]))
    remaining_boxes = boxes.copy()
    solution = []
    
    while remaining_boxes:
        step = []
        available_lifters = set(range(len(lifters)))
        
        # Try to assign boxes in current step
        i = 0
        while i < len(remaining_boxes):
            if not available_lifters:
                break
                
            box_idx, weight = remaining_boxes[i]
            
            # Try to find best lifter combination for current box
            best_combo = find_best_combination(weight, lifter_combinations, available_lifters)
            
            if best_combo:
                step.append((weight, list(best_combo)))
                for lifter in best_combo:
                    available_lifters.remove(lifter)
                remaining_boxes.pop(i)
            else:
                i += 1
        
        if step:
            solution.append(step)
            if len(solution) > 5:
                print("<<<No solution found within 5 steps>>>")
                return
        else:
            # If we can't make progress in a step, no solution exists
            print("<<<No solution found within 5 steps>>>")
            return
    
    # Format solution
    result = []
    for i, step in enumerate(solution, 1):
        result.append(f"Step {i}: {step}")
    print("<<<" + "\n".join(result) + ">>>")

solve_box_lifting()