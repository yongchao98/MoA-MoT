from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_solution():
    boxes = [159, 180, 64, 184, 110, 43, 128, 47, 47, 138, 62, 97]
    lifter_capacities = [81, 98, 62, 46]
    boxes = sorted(list(enumerate(boxes)), key=lambda x: x[1], reverse=True)
    
    solution = []
    remaining_boxes = boxes[:]
    step = 1
    
    while remaining_boxes and step <= 7:
        step_solution = []
        available_lifters = set(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        current_boxes = remaining_boxes[:]
        
        while current_boxes and available_lifters:
            box_idx, box_weight = current_boxes[0]
            assigned = False
            
            # Try combinations of available lifters for current box
            for r in range(1, len(available_lifters) + 1):
                if assigned:
                    break
                for lifter_combo in combinations(available_lifters, r):
                    lifter_sum = sum(lifter_capacities[i] for i in lifter_combo)
                    if lifter_sum >= box_weight:
                        step_solution.append((box_weight, list(lifter_combo)))
                        available_lifters -= set(lifter_combo)
                        current_boxes.pop(0)
                        assigned = True
                        break
            
            if not assigned:
                break
        
        remaining_boxes = current_boxes
        if step_solution:
            solution.append(f"Step {step}: {step_solution}")
        step += 1
    
    if not remaining_boxes:
        print("<<<" + "\n".join(solution) + ">>>")
    else:
        print("No solution found within 7 steps")

find_solution()