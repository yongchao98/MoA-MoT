from itertools import combinations

boxes = [187, 344, 112, 66, 384, 247, 184, 298, 131, 51, 62, 255, 120, 357, 399, 287, 231, 161, 336, 256, 328, 239, 365, 245]
lifters = [165, 85, 52, 82, 114, 75]

# Sort boxes and lifters in descending order
boxes.sort(reverse=True)
lifters.sort(reverse=True)

def can_lift(box, lifters):
    for r in range(1, len(lifters) + 1):
        for combo in combinations(lifters, r):
            if sum(combo) >= box:
                return combo
    return None

def dynamic_programming_lift(boxes, lifters):
    n = len(boxes)
    dp = [float('inf')] * (1 << n)
    dp[0] = 0

    for mask in range(1 << n):
        if dp[mask] < float('inf'):
            available_boxes = [boxes[i] for i in range(n) if not (mask & (1 << i))]
            for box in available_boxes:
                lifter_combo = can_lift(box, lifters)
                if lifter_combo:
                    new_mask = mask
                    for i in range(n):
                        if boxes[i] == box and not (new_mask & (1 << i)):
                            new_mask |= (1 << i)
                            break
                    dp[new_mask] = min(dp[new_mask], dp[mask] + 1)

    return dp[(1 << n) - 1]

min_steps = dynamic_programming_lift(boxes, lifters)

if min_steps <= 11:
    print(f"<<<Solution found in {min_steps} steps>>>")
else:
    print("No solution found within 11 steps.")