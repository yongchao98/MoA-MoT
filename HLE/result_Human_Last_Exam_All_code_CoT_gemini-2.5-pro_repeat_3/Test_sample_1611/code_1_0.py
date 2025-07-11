import numpy as np

def solve_grid_assignments():
    """
    Calculates the number of valid 0/1 assignments for an n x m grid.
    """
    n = 4
    m = 4

    # Step 1 & 2: Generate valid row states for a 1xm grid
    # A valid row has no adjacent 1s.
    valid_rows = []
    for i in range(1 << m):
        # Check if there are any adjacent 1s
        if (i & (i << 1)) == 0:
            valid_rows.append(i)
    
    k = len(valid_rows)
    states = sorted(valid_rows)

    # Step 3: Create the transition matrix T
    T = np.zeros((k, k), dtype=np.int64)
    for i in range(k):
        for j in range(k):
            # Two rows are compatible if they have no 1s in the same column
            if (states[i] & states[j]) == 0:
                T[i, j] = 1

    # Step 4: Calculate assignments for n rows
    # v1 represents the number of ways to form a 1xm grid ending in each state.
    # Initially, each valid state is one way.
    v = np.ones(k, dtype=np.int64)

    # Sequentially calculate v_n for n=2, 3, 4
    v_prev = v
    for i in range(2, n + 2):
        v_curr = v_prev @ T
        if i == n:
            # Step 5: Final calculation for the 4x4 grid
            total_sum = np.sum(v_curr)
            equation_parts = [str(val) for val in v_curr]
            equation_str = " + ".join(equation_parts)
            print(f"{equation_str} = {total_sum}")
            break
        v_prev = v_curr

solve_grid_assignments()