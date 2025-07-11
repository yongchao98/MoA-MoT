import math

def calculate_tb():
    """
    Calculates the maximal Thurston-Bennequin number for the knot
    associated with the given 5x5 grid diagram.
    """
    n = 5
    # O's are at (i, sigma[i]), but since they are at (i,i), sigma is the identity.
    # X's are at (i, tau[i]).
    # We use 1-based indexing for permutations to match the problem statement.
    # A dummy 0 is placed at index 0.
    tau = [0, 4, 5, 1, 2, 3]

    # --- Step 1: Calculate the writhe of the grid diagram G ---
    # tau_inv[j] gives the column i where the X in row j is located.
    tau_inv = [0] * (n + 1)
    for i in range(1, n + 1):
        tau_inv[tau[i]] = i

    writhe = 0
    crossing_locations = []
    # Iterate over all possible grid points (i, j) to find crossings.
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            # A crossing exists at (i, j) if i is strictly between j and tau_inv[j]
            # AND j is strictly between i and tau[i].
            # (assuming sigma is identity).
            cond1 = (j < i < tau_inv[j]) or (tau_inv[j] < i < j)
            cond2 = (i < j < tau[i]) or (tau[i] < j < i)

            if cond1 and cond2:
                # Vertical segment direction: sign(tau[i] - i)
                sign_v = 1 if tau[i] > i else -1
                # Horizontal segment direction: sign(i_of_O - i_of_X) = sign(j - tau_inv[j])
                sign_h = 1 if j > tau_inv[j] else -1
                # Sign of crossing is product of orientations
                crossing_sign = sign_v * sign_h
                writhe += crossing_sign
                crossing_locations.append((i,j))
                
    # --- Step 2: Calculate tb(G) for the given grid ---
    # The formula for the tb number of this specific grid representative is tb(G) = writhe - n
    tb_G = writhe - n
    
    # --- Step 3: Count the number of stabilizations ---
    # A grid is stabilized if there exists a pair i < j such that tau[i] < tau[j]
    # and the rectangle (i, j) x (tau[i], tau[j]) is empty of other markings.
    num_stabilizations = 0
    stabilization_locations = []
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if tau[i] < tau[j]:
                # This pair (i,j) is a potential stabilization. Check if the rectangle is empty.
                is_empty = True
                # Check for O-marks (k,k) in the rectangle's interior.
                for k in range(i + 1, j):
                    if tau[i] < k < tau[j]:
                        is_empty = False
                        break
                if not is_empty:
                    continue

                # Check for X-marks (k, tau[k]) in the rectangle's interior.
                for k in range(1, n + 1):
                    if k != i and k != j: # Exclude the corners of the rectangle
                        if (i < k < j) and (tau[i] < tau[k] < tau[j]):
                            is_empty = False
                            break
                
                if is_empty:
                    num_stabilizations += 1
                    stabilization_locations.append((i,j))

    # --- Step 4: Calculate the maximal tb number ---
    # Each stabilization decreases the tb number by 1. So, TB_max = tb(G) + num_stabilizations.
    tb_max = tb_G + num_stabilizations
    
    # --- Step 5: Output the results ---
    print(f"The grid size is n = {n}.")
    print(f"The permutation tau defined by the X's is: {tau[1:]}")
    print(f"\nStep 1: Calculate the writhe of the grid diagram.")
    print(f"Found {len(crossing_locations)} crossings at locations {crossing_locations}.")
    print(f"The total writhe of the diagram is w(G) = {writhe}.")
    
    print(f"\nStep 2: Calculate the Thurston-Bennequin number for this specific grid, tb(G).")
    print(f"Using the formula tb(G) = w(G) - n:")
    print(f"tb(G) = {writhe} - {n} = {tb_G}")
    
    print(f"\nStep 3: Count stabilizations to find the maximal tb.")
    print(f"Found {num_stabilizations} stabilizations at column pairs {stabilization_locations}.")
    
    print(f"\nStep 4: Calculate the maximal Thurston-Bennequin number, TB(K).")
    print(f"Using the formula TB(K) = tb(G) + (number of stabilizations):")
    print(f"TB(K) = {tb_G} + {num_stabilizations} = {tb_max}")
    print(f"\nThe maximal Thurston-Bennequin number is {tb_max}.")

calculate_tb()