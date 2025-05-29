from heapq import heappush, heappop

def can_lift(box, lifters, used):
    lifter_indices = []
    remaining_weight = box
    for i, lifter in enumerate(lifters):
        if not used[i] and lifter <= remaining_weight:
            lifter_indices.append(i)
            remaining_weight -= lifter
            used[i] = True
        if remaining_weight <= 0:
            return lifter_indices
    for idx in lifter_indices:
        used[idx] = False
    return None

def branch_and_bound(boxes, lifters):
    boxes.sort(reverse=True)
    lifters.sort(reverse=True)
    steps = [[] for _ in range(7)]
    used = [False] * len(lifters)
    
    # Priority queue for branch-and-bound
    pq = []
    heappush(pq, (0, 0, boxes, steps, used))
    
    while pq:
        neg_lifted, step, remaining_boxes, current_steps, current_used = heappop(pq)
        
        if not remaining_boxes:
            return current_steps
        
        if step >= 7:
            continue
        
        for i, box in enumerate(remaining_boxes):
            lifter_indices = can_lift(box, lifters, current_used[:])
            if lifter_indices is not None:
                new_steps = [s[:] for s in current_steps]
                new_steps[step].append((box, lifter_indices))
                new_remaining_boxes = remaining_boxes[:i] + remaining_boxes[i+1:]
                new_used = current_used[:]
                for idx in lifter_indices:
                    new_used[idx] = True
                total_capacity = sum(lifters[idx] for idx in range(len(lifters)) if not new_used[idx])
                heappush(pq, (neg_lifted - 1, step, new_remaining_boxes, new_steps, new_used))
        
        # Move to the next step
        heappush(pq, (neg_lifted, step + 1, remaining_boxes, current_steps, [False] * len(lifters)))
    
    return "Not possible to lift all boxes in 7 steps or less."

boxes = [362, 202, 192, 149, 357, 335, 209, 193, 334, 375, 374, 167, 276, 70, 282, 201, 242, 158, 236, 212]
lifters = [174, 130, 182, 186, 80, 149, 87]

steps = branch_and_bound(boxes, lifters)

output = "<<<"
for i, step in enumerate(steps):
    if step:
        output += f"Step {i+1}: {step}\n"
output += ">>>"

print(output)