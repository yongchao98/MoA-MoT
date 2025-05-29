from itertools import combinations

boxes = [88, 128, 39, 30, 80, 60, 162, 154, 191, 150, 178, 32]
lifters = [95, 111, 115, 35]

# Sort boxes in descending order
boxes.sort(reverse=True)

# Function to check if a combination of lifters can lift a box
def can_lift(box, lifter_indices):
    return sum(lifters[i] for i in lifter_indices) >= box

# Function to find the optimal solution
def find_optimal_solution(boxes, lifters):
    n = len(boxes)
    dp = [None] * (1 << n)
    dp[0] = []

    for mask in range(1 << n):
        if dp[mask] is None:
            continue

        for step in range(5):
            used_lifters = [False] * len(lifters)
            new_mask = mask

            for i in range(n):
                if not (mask & (1 << i)):
                    for r in range(1, len(lifters) + 1):
                        for lifter_indices in combinations(range(len(lifters)), r):
                            if all(not used_lifters[j] for j in lifter_indices) and can_lift(boxes[i], lifter_indices):
                                new_mask |= (1 << i)
                                for j in lifter_indices:
                                    used_lifters[j] = True
                                break

            if new_mask != mask:
                if dp[new_mask] is None or len(dp[new_mask]) > len(dp[mask]) + 1:
                    dp[new_mask] = dp[mask] + [step]

    return dp[(1 << n) - 1]

# Find the solution
solution = find_optimal_solution(boxes, lifters)

# Print the solution
if solution and len(solution) <= 5:
    output = []
    for i, step in enumerate(solution):
        output.append(f"Step {i + 1}: {step}")
    print("<<<" + "\n".join(output) + ">>>")
else:
    print("No solution found within 5 steps.")