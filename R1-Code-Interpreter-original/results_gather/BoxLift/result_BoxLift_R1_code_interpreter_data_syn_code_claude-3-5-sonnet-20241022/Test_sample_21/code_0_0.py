from itertools import combinations
import copy

def can_lift(lifters, weight):
    return sum(lifters) >= weight

def find_solution():
    boxes = [38, 97, 21, 92, 195, 187, 95, 140, 37, 85, 25, 165]
    lifter_capacities = [96, 32, 88, 48]
    boxes.sort(reverse=True)  # Sort boxes in descending order
    
    remaining_boxes = boxes.copy()
    steps = []
    
    while remaining_boxes and len(steps) < 6:
        step = []
        available_lifters = list(range(len(lifter_capacities)))
        
        # Try to assign boxes to lifters
        while available_lifters and remaining_boxes:
            assigned = False
            
            # Try different combinations of lifters for current heaviest box
            for box_idx, box_weight in enumerate(remaining_boxes):
                for r in range(1, len(available_lifters) + 1):
                    for lifter_combo in combinations(available_lifters, r):
                        total_capacity = sum(lifter_capacities[i] for i in lifter_combo)
                        
                        if total_capacity >= box_weight:
                            step.append((box_weight, list(lifter_combo)))
                            remaining_boxes.pop(box_idx)
                            for lifter in lifter_combo:
                                available_lifters.remove(lifter)
                            assigned = True
                            break
                    if assigned:
                        break
                if assigned:
                    break
            
            if not assigned:
                break
        
        if step:
            steps.append(step)
        else:
            break
    
    # Format and print the solution
    if not remaining_boxes and len(steps) <= 6:
        result = ""
        for i, step in enumerate(steps, 1):
            result += f"Step {i}: {step}\n"
        print(result.strip())
    else:
        print("No solution found within 6 steps")

find_solution()