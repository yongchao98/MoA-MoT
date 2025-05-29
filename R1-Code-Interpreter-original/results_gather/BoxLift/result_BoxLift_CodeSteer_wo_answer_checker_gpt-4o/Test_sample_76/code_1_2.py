from functools import lru_cache

boxes = [94, 177, 88, 131, 87, 113, 148, 86, 141, 83, 86, 73]
lifters = [55, 92, 35, 48, 96]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

@lru_cache(None)
def can_lift(box, available_lifters):
    # Try to find a combination of lifters to lift the box
    if box <= 0:
        return []
    
    for i, lifter in enumerate(available_lifters):
        if lifter <= box:
            remaining_lifters = available_lifters[:i] + available_lifters[i+1:]
            result = can_lift(box - lifter, remaining_lifters)
            if result is not None:
                return [i] + result
    
    return None

steps = []
while boxes:
    step = []
    available_lifters = tuple(lifters)
    
    for box in boxes[:]:
        lifter_indices = can_lift(box, available_lifters)
        if lifter_indices is not None:
            step.append((box, lifter_indices))
            for index in lifter_indices:
                available_lifters = available_lifters[:index] + available_lifters[index+1:]
            boxes.remove(box)
    
    steps.append(step)

# Print the steps
output = []
for i, step in enumerate(steps):
    output.append(f"Step {i+1}: {step}")

print("<<<" + "\n".join(output) + ">>>")