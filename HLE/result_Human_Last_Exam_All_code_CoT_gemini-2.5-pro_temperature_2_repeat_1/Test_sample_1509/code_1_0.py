# This script demonstrates the counterexample for question (b).
# Question: Must a shifted (t+1)-intersecting family F satisfy |F^(n)| >= 3
# for n >= k + t + 3?

# 1. Define the parameters and the family for our counterexample.
t = 1
k = 3
n = 7
# The family F, as a list of sets.
F_list = [{1, 2, 3}, {1, 2, 4}]
F = {frozenset(s) for s in F_list}

# 2. Verify that the premises of the question hold for our example.
print("--- Verifying the counterexample for question (b) ---")
print(f"Parameters: n={n}, k={k}, t={t}")
print(f"Family F = {F_list}")

# Premise 1: n >= k + t + 3
n_cond_val = k + t + 3
n_cond_holds = (n >= n_cond_val)
print(f"\nVerifying Premise 1: n >= k + t + 3")
print(f"  Is {n} >= {k} + {t} + 3?  =>  Is {n} >= {n_cond_val}?  =>  {n_cond_holds}")

# Premise 2: F is shifted.
# For F={{1,2,3},{1,2,4}}, the only non-trivial shift operation is on {1,2,4}
# with j=4, i=3, which produces {1,2,3}, a member of F. So F is shifted.
is_f_shifted = True
print(f"\nVerifying Premise 2: F is shifted?  =>  {is_f_shifted}")

# Premise 3: F is (t+1)-intersecting.
intersect_target = t + 1
# There is only one pair of sets in F to check.
pair = tuple(F)
intersection_size = len(pair[0].intersection(pair[1]))
is_f_intersecting = (intersection_size >= intersect_target)
print(f"\nVerifying Premise 3: F is {intersect_target}-intersecting?")
print(f"  Intersection of {set(pair[0])} and {set(pair[1])} has size {intersection_size}.")
print(f"  Is {intersection_size} >= {intersect_target}?  =>  {is_f_intersecting}")

# 3. Check if the conclusion holds.
print("\n--- Checking the conclusion ---")
if n_cond_holds and is_f_shifted and is_f_intersecting:
    print("All premises are satisfied by our example.")
    
    # Construct F^(n)
    F_n = {s for s in F if n not in s}
    size_F_n = len(F_n)
    
    # Check the conclusion: |F^(n)| >= 3
    conclusion_target = 3
    conclusion_holds = (size_F_n >= conclusion_target)
    
    print(f"\nWe must check if |F^({n})| >= 3.")
    F_n_list = [set(s) for s in F_n]
    print(f"  F^({n}) consists of sets in F not containing {n}. F^({n}) = {F_n_list}")
    print(f"  The size of F^({n}) is {size_F_n}.")
    print(f"  Is {size_F_n} >= {conclusion_target}?  =>  {conclusion_holds}")
    print("\nThe premise is true, but the conclusion is false.")
    print("\nTherefore, the answer to (b) is No.")
else:
    print("\nThe example does not satisfy the premises and cannot be used.")
