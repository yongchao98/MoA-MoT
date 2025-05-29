import math

def find_feasible_input():
    x1, y1 = 39, 73
    # Start with some initial guesses for x2 and y2
    x2_candidates = [38, 37, 36, 35]
    y2_candidates = [73, 74, 75, 76]
    
    for x2 in x2_candidates:
        for y2 in y2_candidates:
            res1 = y1 * math.log(x1)
            res2 = y2 * math.log(x2)
            if res1 > res2:
                return {"x1": x1, "y1": y1, "x2": x2, "y2": y2}

# Execute the function and print the result
print(find_feasible_input())