boxes = [50, 219, 120, 245, 204, 172, 274, 225, 94, 270, 158, 253, 174, 40, 252, 124]
lifters = [54, 61, 106, 60, 63]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box, lifters, used):
    total_capacity = 0
    lifter_indices = []
    for i, capacity in enumerate(lifters):
        if not used[i] and total_capacity + capacity <= box:
            total_capacity += capacity
            lifter_indices.append(i)
            if total_capacity == box:
                break
    return total_capacity >= box, lifter_indices

def branch_and_bound(boxes, lifters, step, max_steps, used):
    if not boxes:
        return step
    if step >= max_steps:
        return float('inf')
    
    state = (tuple(boxes), tuple(used))
    if state in memo:
        return memo[state]
    
    min_steps = float('inf')
    for i, box in enumerate(boxes):
        can_lift_box, lifter_indices = can_lift(box, lifters, used)
        if can_lift_box:
            new_used = used[:]
            for idx in lifter_indices:
                new_used[idx] = True
            remaining_boxes = boxes[:i] + boxes[i+1:]
            result = branch_and_bound(remaining_boxes, lifters, step + 1, max_steps, new_used)
            min_steps = min(min_steps, result)
    
    memo[state] = min_steps
    return min_steps

# Initialize used lifters
used_lifters = [False] * len(lifters)
min_steps = branch_and_bound(boxes, lifters, 0, 12, used_lifters)

print(f"Minimum steps required: {min_steps}")