from functools import lru_cache

boxes = [299, 338, 270, 207, 76, 195, 101, 89, 277, 236, 240, 342, 47, 207, 220, 349, 185, 124, 309, 112, 266, 278, 162, 119]
lifters = [198, 55, 76, 137, 163, 117, 146]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization cache
@lru_cache(None)
def can_lift(remaining_boxes, available_lifters):
    if not remaining_boxes:
        return True, []

    box = remaining_boxes[0]
    for i in range(len(available_lifters)):
        if available_lifters[i] >= box:
            new_lifters = available_lifters[:i] + available_lifters[i+1:]
            can_lift_rest, steps = can_lift(tuple(remaining_boxes[1:]), tuple(new_lifters))
            if can_lift_rest:
                return True, [(box, [i])] + steps

        # Try combining lifters
        for j in range(i+1, len(available_lifters)):
            if available_lifters[i] + available_lifters[j] >= box:
                new_lifters = available_lifters[:i] + available_lifters[i+1:j] + available_lifters[j+1:]
                can_lift_rest, steps = can_lift(tuple(remaining_boxes[1:]), tuple(new_lifters))
                if can_lift_rest:
                    return True, [(box, [i, j])] + steps

    return False, []

# Find the solution
success, steps = can_lift(tuple(boxes), tuple(lifters))

# Output the steps
if success:
    output = []
    for i, step in enumerate(steps):
        output.append(f"Step {i+1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within the constraints.")