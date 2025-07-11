import sys

def solve():
    """
    Solves for the maximal Thurston-Bennequin number of the knot from the grid diagram.
    """
    # Step 1: Define the grid parameters from the problem description.
    n = 5
    o_positions = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    x_positions = [(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)]
    
    # Extract the permutation pi from the X positions.
    # The permutation maps column i to the row of the X in that column.
    pi = [0] * n # Use 0-based indexing for the list, for pi[i-1]
    for c, r in x_positions:
        pi[c-1] = r

    # Step 2: Identify the knot type by finding the parameter q for the torus knot T(n, q).
    # The permutation for T(n, q) is pi_q(i) = ((i-1 + q) % n) + 1 for i in {1,...,n}.
    p = n
    found_q = -1
    for q_candidate in range(1, n):
        is_match = True
        for i in range(1, n + 1):
            expected_pi_i = ((i - 1 + q_candidate) % n) + 1
            if pi[i-1] != expected_pi_i:
                is_match = False
                break
        if is_match:
            found_q = q_candidate
            break

    if found_q == -1:
        print("The knot is not a simple torus knot of the form T(n,q).")
        return

    q = found_q
    
    # Step 3: Use the formula for the maximal Thurston-Bennequin number.
    # tb_max(T(p,q)) = p*q - p - q
    tb_max = p * q - p - q

    # Step 4: Output the result and the equation.
    print(f"The grid diagram corresponds to the torus knot T({p},{q}).")
    print("The formula for the maximal Thurston-Bennequin number for a torus knot T(p,q) is: pq - p - q.")
    print(f"For T({p},{q}), the calculation is: {p} * {q} - {p} - {q} = {tb_max}")
    print(f"The maximal Thurston-Bennequin number is {tb_max}.")

solve()