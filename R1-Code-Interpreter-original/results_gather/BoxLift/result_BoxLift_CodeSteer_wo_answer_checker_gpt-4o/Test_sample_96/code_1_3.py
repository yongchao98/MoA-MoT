boxes = [169, 106, 153, 20, 138, 80, 94, 190, 142, 196, 186, 87]
lifters = [55, 85, 43, 100, 47]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

from functools import lru_cache

@lru_cache(None)
def min_steps(remaining_boxes, available_lifters):
    if not remaining_boxes:
        return 0
    
    min_steps_needed = float('inf')
    for i, box in enumerate(remaining_boxes):
        for j in range(1 << len(available_lifters)):
            total_capacity = 0
            lifter_indices = []
            for k in range(len(available_lifters)):
                if j & (1 << k):
                    total_capacity += available_lifters[k]
                    lifter_indices.append(k)
            
            if total_capacity >= box:
                new_remaining_boxes = remaining_boxes[:i] + remaining_boxes[i+1:]
                new_available_lifters = tuple(l for idx, l in enumerate(available_lifters) if idx not in lifter_indices)
                steps = 1 + min_steps(new_remaining_boxes, new_available_lifters)
                min_steps_needed = min(min_steps_needed, steps)
    
    return min_steps_needed

# Convert lists to tuples for hashing in lru_cache
remaining_boxes = tuple(boxes)
available_lifters = tuple(lifters)

# Calculate the minimum steps
steps_needed = min_steps(remaining_boxes, available_lifters)

print(f"Minimum steps needed: {steps_needed}")