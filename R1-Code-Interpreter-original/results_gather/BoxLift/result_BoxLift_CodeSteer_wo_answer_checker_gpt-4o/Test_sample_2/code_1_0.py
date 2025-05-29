import heapq

boxes = [66, 299, 90, 133, 52, 245, 57, 62, 71, 217, 117, 193, 204, 84, 224, 51]
lifters = [123, 155, 65, 92, 92]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Create a priority queue for lifters with negative capacities for max-heap behavior
lifter_queue = [(-capacity, i) for i, capacity in enumerate(lifters)]
heapq.heapify(lifter_queue)

def find_combination(box, lifter_queue, used):
    # Try to find a single lifter or a combination of lifters to lift the box
    best_combination = None
    best_capacity_used = 0
    available_lifters = []

    while lifter_queue:
        capacity, index = heapq.heappop(lifter_queue)
        capacity = -capacity
        if used[index]:
            available_lifters.append((-capacity, index))
            continue
        if capacity >= box:
            available_lifters.append((-capacity, index))
            return [index]
        for _, (next_capacity, next_index) in enumerate(available_lifters):
            combined_capacity = capacity + (-next_capacity)
            if combined_capacity >= box and combined_capacity > best_capacity_used:
                best_combination = [index, next_index]
                best_capacity_used = combined_capacity
        available_lifters.append((-capacity, index))

    # Push back all available lifters to the queue
    for lifter in available_lifters:
        heapq.heappush(lifter_queue, lifter)

    return best_combination

def backtrack(boxes, lifter_queue, step, steps, used):
    if not boxes:
        return True
    if step >= 6:
        return False

    current_step = []
    remaining_boxes = boxes[:]
    for box in boxes:
        lifter_indices = find_combination(box, lifter_queue, used)
        if lifter_indices is not None:
            current_step.append((box, lifter_indices))
            for index in lifter_indices:
                used[index] = True
            remaining_boxes.remove(box)

    if current_step:
        steps.append(current_step)
        if backtrack(remaining_boxes, lifter_queue, step + 1, steps, [False] * len(lifters)):
            return True
        steps.pop()

    return False

steps = []
used = [False] * len(lifters)
if backtrack(boxes, lifter_queue, 0, steps, used):
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 6 steps.")