import numpy as np

def solve():
    """
    Calculates the number of valid 0/1 assignments for a 4x4 grid
    based on the given implication rules.
    """
    n = 4  # grid height
    m = 4  # grid width

    # Step 1: Find all valid states for a single column.
    # A state is valid if no two vertically adjacent cells are 1.
    valid_column_states = []
    for i in range(1 << n):
        # Check for adjacent 1s in the binary representation of i
        if not any((i >> j) & 1 and (i >> (j + 1)) & 1 for j in range(n - 1)):
            valid_column_states.append(i)

    num_states = len(valid_column_states)
    state_map = {state: i for i, state in enumerate(valid_column_states)}

    # Step 2: Build the transition matrix T.
    # T[i, j] = 1 if column state i can be adjacent to column state j.
    # This is true if they have no 1s at the same vertical position.
    T = np.zeros((num_states, num_states), dtype=object)
    for i in range(num_states):
        for j in range(num_states):
            state_i = valid_column_states[i]
            state_j = valid_column_states[j]
            if (state_i & state_j) == 0:
                T[i, j] = 1

    # Step 3: Use the transfer matrix method to count assignments.
    # u_k[i] = number of valid n x k grids ending with state i.
    # u_1 is all ones, as any valid state can be the first column.
    u = np.ones(num_states, dtype=object)

    # We need to compute T^(m-1) * u.
    # Since T is symmetric, u_k = T @ u_{k-1}.
    T_power = np.linalg.matrix_power(T, m - 1)
    final_counts = T_power.dot(u)

    # The total number of assignments is the sum of all ways.
    total_assignments = sum(final_counts)
    
    # Step 4: Output the final calculation.
    # The problem asks to output each number in the final equation.
    # This is the sum of the elements in the final vector.
    sum_str = " + ".join(map(str, final_counts))
    print(f"The number of valid assignments for each possible last column state are: {final_counts}")
    print("\nThe final equation for the total number of assignments is:")
    print(f"{sum_str} = {total_assignments}")


solve()