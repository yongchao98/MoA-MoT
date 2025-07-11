import numpy as np

def solve_knot_properties():
    """
    Analyzes a 5x5 grid diagram to find the maximal Thurston-Bennequin number
    of its associated knot.
    """
    n = 5
    # O's are at (i,i), so pi_O is the identity.
    # X's are at (i, pi_X(i)).
    # We are given X's at positions (1,4), (2,5), (3,1), (4,2), (5,3).
    # Python uses 0-based indexing, so we convert.
    # pi_X maps 0 to 3, 1 to 4, 2 to 0, 3 to 1, 4 to 2.
    pi = {i + 1: val for i, val in enumerate([4, 5, 1, 2, 3])}
    
    # Calculate the inverse permutation pi_inv
    pi_inv = {pi[i]: i for i in pi}
    
    print("Step 1: Define permutations from grid coordinates.")
    print(f"pi_X = {pi}")
    print(f"pi_X_inv = {pi_inv}\n")
    
    writhe = 0
    crossings = []
    
    # Step 2: Calculate the writhe by finding all crossings and their signs.
    # A crossing exists at (i, j) if the vertical segment for column i and the
    # horizontal segment for row j intersect.
    print("Step 2: Calculate the writhe of the knot diagram.")
    print("Finding all crossings and their signs...")
    
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if i == j:
                continue
            
            # Check if j is on the vertical path for column i
            is_on_vert_path = (pi[i] > i and i < j < pi[i]) or \
                              (pi[i] < i and pi[i] < j < i)
                              
            # Check if i is on the horizontal path for row j
            is_on_horz_path = (pi_inv[j] > j and j < i < pi_inv[j]) or \
                              (pi_inv[j] < j and pi_inv[j] < i < j)
            
            if is_on_vert_path and is_on_horz_path:
                # Orientation of vertical segment: sign(pi[i] - i)
                sign_v = np.sign(pi[i] - i)
                # Orientation of horizontal segment: sign(j - pi_inv[j])
                sign_h = np.sign(j - pi_inv[j])
                # Crossing sign is the product
                sign = sign_v * sign_h
                crossings.append({'pos': (i, j), 'sign': int(sign)})
                writhe += sign

    print("Crossings found at (column, row) with sign:")
    equation_parts = []
    for crossing in crossings:
        sign_str = "+1" if crossing['sign'] > 0 else "-1"
        print(f"  {crossing['pos']}: {sign_str}")
        equation_parts.append(f"({sign_str})")

    equation = " + ".join(equation_parts)
    print(f"\nThe total writhe is the sum of these signs:")
    print(f"Writhe = {equation} = {writhe}\n")

    # Step 3: Find the number of components by cycle decomposition of pi
    print("Step 3: Determine the number of components of the link.")
    num_components = 0
    visited = set()
    cycles = []
    for i in range(1, n + 1):
        if i not in visited:
            num_components += 1
            cycle = []
            curr = i
            while curr not in visited:
                visited.add(curr)
                cycle.append(curr)
                curr = pi[curr]
            cycles.append(cycle)

    print(f"The permutation pi decomposes into {num_components} cycle(s):")
    for cycle in cycles:
        print(f"  {tuple(cycle)}")
    
    print("\nSince there is only 1 component, the diagram represents a knot.\n")

    # Step 4: Identify the knot and find its maximal TB number.
    print("Step 4: Identify the knot and find its maximal Thurston-Bennequin number.")
    print(f"The diagram has a crossing number of {len(crossings)} and a writhe of {writhe}.")
    print("A knot with 3 crossings and writhe +3 is the right-handed trefoil knot (3_1).")
    print("The maximal Thurston-Bennequin number for the right-handed trefoil knot is a known invariant.")
    max_tb = 1
    print(f"TB(3_1) = {max_tb}")
    
    return max_tb

if __name__ == '__main__':
    final_answer = solve_knot_properties()
    print(f"\nThe final answer is {final_answer}")
    print(f'<<<{final_answer}>>>')
