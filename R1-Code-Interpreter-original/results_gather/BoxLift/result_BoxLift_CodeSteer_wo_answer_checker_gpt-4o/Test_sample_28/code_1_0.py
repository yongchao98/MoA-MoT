from itertools import combinations
from queue import PriorityQueue

def can_lift(box, lifters, used):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(range(len(lifters)), r):
            if all(not used[i] for i in combo) and sum(lifters[i] for i in combo) >= box:
                return combo
    return None

def branch_and_bound(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    pq = PriorityQueue()
    pq.put((0, boxes, []))  # (step, remaining_boxes, steps_taken)
    
    while not pq.empty():
        step, remaining_boxes, steps_taken = pq.get()
        
        if not remaining_boxes:
            return steps_taken
        
        if step >= 7:
            continue
        
        used = [False] * len(lifters)
        current_step = []
        new_remaining_boxes = remaining_boxes[:]
        
        for box in remaining_boxes:
            lifter_indices = can_lift(box, lifters, used)
            if lifter_indices is not None:
                current_step.append((box, list(lifter_indices)))
                new_remaining_boxes.remove(box)
                for i in lifter_indices:
                    used[i] = True
        
        if current_step:
            new_steps_taken = steps_taken + [current_step]
            pq.put((step + 1, new_remaining_boxes, new_steps_taken))
    
    return "Cannot lift all boxes in 7 steps or less."

boxes = [71, 72, 167, 64, 339, 45, 243, 53, 52, 280, 164, 322, 362, 75, 158, 142, 150, 353, 182, 282]
lifters = [139, 55, 172, 82, 129, 59, 130]

steps = branch_and_bound(boxes, lifters)

output = "<<<"
if isinstance(steps, str):
    output += steps
else:
    for i, step in enumerate(steps):
        output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)