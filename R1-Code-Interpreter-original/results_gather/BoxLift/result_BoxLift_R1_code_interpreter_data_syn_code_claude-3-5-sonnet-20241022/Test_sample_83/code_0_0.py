import itertools
from typing import List, Tuple, Dict

def can_lift(lifters: List[int], weight: int) -> List[List[int]]:
    valid_combinations = []
    for r in range(1, len(lifters) + 1):
        for combo in itertools.combinations(range(len(lifters)), r):
            if sum(lifters[i] for i in combo) >= weight:
                valid_combinations.append(list(combo))
    return valid_combinations

def solve_box_lifting(boxes: List[int], lifters: List[int], max_steps: int = 7) -> List[List[Tuple[int, List[int]]]]:
    boxes = sorted(boxes, reverse=True)  # Sort boxes by weight (heaviest first)
    n_boxes = len(boxes)
    best_solution = None
    min_steps = float('inf')
    
    def try_assignment(box_idx: int, current_step: int, steps: List[List[Tuple[int, List[int]]]], 
                      available_lifters: List[bool]) -> None:
        nonlocal best_solution, min_steps
        
        if current_step >= min_steps or current_step > max_steps:
            return
            
        if box_idx >= n_boxes:
            if current_step < min_steps:
                min_steps = current_step
                best_solution = [step[:] for step in steps]
            return
            
        # Try to lift current box
        weight = boxes[box_idx]
        available_lifter_indices = [i for i, available in enumerate(available_lifters) if available]
        valid_combinations = can_lift([lifters[i] for i in available_lifter_indices], weight)
        
        for combo in valid_combinations:
            actual_lifters = [available_lifter_indices[i] for i in combo]
            
            # Try adding to current step
            if any(available_lifters):
                temp_available = available_lifters[:]
                for lifter_idx in actual_lifters:
                    temp_available[lifter_idx] = False
                
                if len(steps) <= current_step:
                    steps.append([])
                steps[current_step].append((weight, actual_lifters))
                
                try_assignment(box_idx + 1, current_step, steps, temp_available)
                
                steps[current_step].pop()
                if not steps[current_step]:
                    steps.pop()
            
            # Try starting new step
            if len(steps) <= current_step + 1:
                steps.append([])
            steps[current_step + 1].append((weight, actual_lifters))
            new_available = [True] * len(lifters)
            for lifter_idx in actual_lifters:
                new_available[lifter_idx] = False
            
            try_assignment(box_idx + 1, current_step + 1, steps, new_available)
            
            steps[current_step + 1].pop()
            if not steps[current_step + 1]:
                steps.pop()

    initial_steps = []
    initial_available = [True] * len(lifters)
    try_assignment(0, -1, initial_steps, initial_available)
    
    return best_solution

# Input data
boxes = [280, 392, 142, 400, 213, 135, 295, 142, 55, 371, 151, 50, 42, 281, 283, 46, 383, 352, 204, 111]
lifters = [143, 62, 162, 128, 141, 159]

# Solve and format output
solution = solve_box_lifting(boxes, lifters)
output = []
for step_idx, step in enumerate(solution, 1):
    output.append(f"Step {step_idx}: {step}")

print("<<<" + "\n".join(output) + ">>>")