import math

def solve_hypersphere_problem():
    """
    Solves the hypersphere point placement problem based on known geometric bounds.
    """
    # n: number of points
    n = 15
    # d: number of dimensions of the hypersphere
    d = 8

    # The problem is to find K, the minimized value of the maximum number of points
    # that can appear in any single closed hyper-hemisphere.
    # K = min_{arrangement P} max_{hemisphere H} |P intersect H|

    # --- Step 1: Calculate the lower bound for K ---
    # For any arrangement of n points, and any hyperplane through the origin,
    # at least one of the two resulting closed hemispheres must contain at least
    # ceil(n/2) points. Therefore, the minimized maximum K cannot be less
    # than this value.
    lower_bound = math.ceil(n / 2)

    # --- Step 2: Calculate the upper bound for K ---
    # A result in discrete geometry (Wu's inequality) states that it is always
    # possible to find an arrangement of n points in d dimensions such that
    # no closed hemisphere contains more than n - d + 1 points (for n >= d+1).
    # This provides an upper bound on the value of K.
    upper_bound = n - d + 1

    # --- Step 3: Combine the bounds to find the solution ---
    print(f"The problem is to place n={n} points in a d={d} dimensional hypersphere.")
    print("The goal is to minimize the maximum number of points in any closed hemisphere.\n")

    print("Step 1: Establishing a lower bound.")
    print(f"For any arrangement, there is always a hemisphere with at least ceil(n/2) points.")
    print(f"Lower Bound = ceil({n}/2) = {int(lower_bound)}\n")

    print("Step 2: Establishing an upper bound.")
    print(f"There always exists an arrangement where no hemisphere has more than n - d + 1 points.")
    print(f"Upper Bound = {n} - {d} + 1 = {upper_bound}\n")
    
    print("Step 3: Conclusion.")
    if lower_bound == upper_bound:
        result = upper_bound
        print("Since the theoretical lower bound and upper bound are equal, we have found the exact value.")
        print("The minimized maximum number of points is given by the equation:")
        print(f"{n} - {d} + 1 = {result}")
    else:
        print(f"The solution lies between {int(lower_bound)} and {upper_bound}.")

solve_hypersphere_problem()