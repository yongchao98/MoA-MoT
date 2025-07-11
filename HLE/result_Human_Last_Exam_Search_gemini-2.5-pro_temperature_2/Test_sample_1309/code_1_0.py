import math

def combinations(n, k):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def R(n, d):
    sum_val = 0
    for k in range(d + 1):
        sum_val += combinations(n, k)
    return sum_val

target = 538902664255516

for d in range(2, 60):
    # Estimate n
    try:
        n_est = (target * math.factorial(d))**(1/d)
    except OverflowError:
        # factorial can become very large
        continue

    n_low = int(n_est) - 5
    n_high = int(n_est) + 5

    for n_test in range(n_low, n_high):
        if n_test < d:
            continue
        val = R(n_test, d)
        if val == target:
            print(f"Found a match: d = {d}, n = {n_test}")
            break
        elif val > target and R(n_test - 1, d) < target:
            # We are close, maybe there is no integer solution for this d
            # no need to continue for this d if we already passed the target
            break
