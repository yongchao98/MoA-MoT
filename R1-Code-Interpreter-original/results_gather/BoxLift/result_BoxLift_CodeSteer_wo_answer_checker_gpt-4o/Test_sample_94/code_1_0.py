import heapq

boxes = [187, 344, 112, 66, 384, 247, 184, 298, 131, 51, 62, 255, 120, 357, 399, 287, 231, 161, 336, 256, 328, 239, 365, 245]
lifters = [165, 85, 52, 82, 114, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def greedy_initial_solution(boxes, lifters):
    steps = 0
    remaining_boxes = boxes[:]
    while remaining_boxes:
        used_lifters = [False] * len(lifters)
        for box in remaining_boxes[:]:
            total_capacity = 0
            for i, lifter in enumerate(lifters):
                if not used_lifters[i] and total_capacity < box:
                    total_capacity += lifter
                    used_lifters[i] = True
                if total_capacity >= box:
                    remaining_boxes.remove(box)
                    break
        steps += 1
    return steps

def branch_and_bound(boxes, lifters):
    initial_steps = greedy_initial_solution(boxes, lifters)
    best_steps = initial_steps
    pq = [(0, 0, boxes)]  # (current_steps, current_mask, remaining_boxes)

    while pq:
        current_steps, _, remaining_boxes = heapq.heappop(pq)
        if current_steps >= best_steps:
            continue
        if not remaining_boxes:
            best_steps = min(best_steps, current_steps)
            continue

        used_lifters = [False] * len(lifters)
        for box in remaining_boxes[:]:
            total_capacity = 0
            for i, lifter in enumerate(lifters):
                if not used_lifters[i] and total_capacity < box:
                    total_capacity += lifter
                    used_lifters[i] = True
                if total_capacity >= box:
                    new_remaining_boxes = remaining_boxes[:]
                    new_remaining_boxes.remove(box)
                    heapq.heappush(pq, (current_steps + 1, 0, new_remaining_boxes))
                    break

    return best_steps

min_steps = branch_and_bound(boxes, lifters)

if min_steps <= 11:
    print(f"<<<Solution found in {min_steps} steps>>>")
else:
    print("No solution found within 11 steps.")