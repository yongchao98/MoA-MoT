import sys

def solve_ck(a_list):
    """
    Calculates c_k for a given list of coefficients a_1, ..., a_k.
    The formula is c_k = p_k * q_{k-1}, where p_k/q_k are convergents of [a_1, ..., a_k].
    """
    k = len(a_list)
    if k < 2:
        print("Error: k must be >= 2.", file=sys.stderr)
        return

    # Use 1-based indexing for a for clarity, matching the problem statement
    a = [0] + a_list 

    # Calculate p_k = N[a_1, ..., a_k]
    # Recurrence: p_n = a_n * p_{n-1} + p_{n-2}
    # Initial values: p_0 = 1, p_1 = a_1
    p = [0] * (k + 1)
    if k >= 0:
      p[0] = 1
    if k >= 1:
      p[1] = a[1]
    for i in range(2, k + 1):
        p[i] = a[i] * p[i - 1] + p[i - 2]
    p_k = p[k]

    # Calculate q_{k-1} = N[a_2, ..., a_{k-1}]
    # This is the (k-1)-th denominator of [a_1, ..., a_k]
    # Recurrence: q_n = a_n * q_{n-1} + q_{n-2}
    # Initial values: q_0 = 0, q_1 = 1
    q = [0] * (k + 1)
    if k >= 0:
      q[0] = 0
    if k >= 1:
      q[1] = 1
    for i in range(2, k): # We only need to go up to k-1
        q[i] = a[i] * q[i - 1] + q[i - 2]
    q_k_minus_1 = q[k - 1]

    # Calculate c_k
    c_k = p_k * q_k_minus_1
    
    # Output the result
    print(f"For k = {k} and a = {a_list}:")
    print(f"p_{k} = N{a_list} = {p_k}")
    print(f"q_{k-1} = N{a_list[1:-1]} = {q_k_minus_1}")
    print(f"c_{k} = p_{k} * q_{k-1} = {p_k} * {q_k_minus_1} = {c_k}")

# Example from the thought process: k=3, a=(2, 3, 4)
print("--- Example 1 ---")
solve_ck([2, 3, 4])

# Example from the thought process: k=2, a=(a1, a2)
# Let a1=5, a2=3
print("\n--- Example 2 ---")
solve_ck([5, 3])

# Another example: k=4, a=(1, 2, 3, 4)
print("\n--- Example 3 ---")
solve_ck([1, 2, 3, 4])
