from itertools import combinations
import copy

def can_lift_box(lifters_combo, box_weight, lifter_capacities):
    return sum(lifter_capacities[i] for i in lifters_combo) >= box_weight

def find_solution():
    boxes = [58, 45, 94, 26, 49, 153, 61, 23, 170, 143, 75, 166]
    lifter_capacities = [57, 61, 104, 98, 70]
    boxes.sort(reverse=True)  # Sort boxes in descending order
    
    solution = []
    remaining_boxes = boxes[:]
    
    while remaining_boxes and len(solution) < 4:
        step = []
        available_lifters = set(range(len(lifter_capacities)))
        current_boxes = remaining_boxes[:]
        
        while current_boxes and available_lifters:
            box = current_boxes[0]
            found_lifters = None
            
            # Try combinations of available lifters
            for r in range(1, len(available_lifters) + 1):
                if found_lifters:
                    break
                for lifters in combinations(available_lifters, r):
                    if can_lift_box(lifters, box, lifter_capacities):
                        found_lifters = list(lifters)
                        break
            
            if found_lifters:
                step.append((box, found_lifters))
                available_lifters -= set(found_lifters)
                current_boxes.remove(box)
                remaining_boxes.remove(box)
        
        if step:
            solution.append(step)
    
    # Format the solution
    if not remaining_boxes and len(solution) <= 4:
        result = ""
        for i, step in enumerate(solution, 1):
            result += f"Step {i}: {step}\n"
        return result.strip()
    return "No solution found"

print(find_solution())