boxes = [55, 77, 274, 297, 70, 312, 108, 84, 340, 273, 157, 386, 95, 223, 311, 306, 287, 91, 220, 296, 320, 237, 97, 73]
lifters = [167, 138, 133, 87, 142, 167, 140]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

# Memoization dictionary
memo = {}

def can_lift(box_index, used_lifters, steps):
    if box_index == len(boxes):
        return steps <= 7
    if steps > 7:
        return False
    state = (box_index, tuple(used_lifters))
    if state in memo:
        return memo[state]
    
    box_weight = boxes[box_index]
    for i in range(len(lifters)):
        if not used_lifters[i]:
            used_lifters[i] = True
            if lifters[i] >= box_weight:
                if can_lift(box_index + 1, used_lifters, steps + 1):
                    memo[state] = True
                    return True
            else:
                for j in range(i + 1, len(lifters)):
                    if not used_lifters[j] and lifters[i] + lifters[j] >= box_weight:
                        used_lifters[j] = True
                        if can_lift(box_index + 1, used_lifters, steps + 1):
                            memo[state] = True
                            return True
                        used_lifters[j] = False
            used_lifters[i] = False
    memo[state] = False
    return False

def find_solution():
    used_lifters = [False] * len(lifters)
    if can_lift(0, used_lifters, 0):
        return "Solution found within 7 steps."
    else:
        return "It's not possible to lift all boxes in 7 or fewer steps."

print(find_solution())