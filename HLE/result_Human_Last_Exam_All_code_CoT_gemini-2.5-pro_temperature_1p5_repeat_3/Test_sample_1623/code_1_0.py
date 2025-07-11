import numpy as np

def solve_thurston_bennequin():
    """
    Calculates the Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # Step 0: Define grid parameters from the problem description
    n = 5
    # O positions (1-based): [(1,1), (2,2), (3,3), (4,4), (5,5)]
    # X positions (1-based): [(1,4), (2,5), (3,1), (4,2), (5,3)]

    # Convert to 0-based for list indexing
    o_pos_0based = [(i - 1, j - 1) for i, j in [(1,1), (2,2), (3,3), (4,4), (5,5)]]
    x_pos_0based = [(i - 1, j - 1) for i, j in [(1,4), (2,5), (3,1), (4,2), (5,3)]]

    # Step 1: Create permutation maps
    pi_O = np.zeros(n, dtype=int)
    pi_X = np.zeros(n, dtype=int)
    sigma_O = np.zeros(n, dtype=int)
    sigma_X = np.zeros(n, dtype=int)

    for i, j in o_pos_0based:
        pi_O[i] = j
        sigma_O[j] = i
    for i, j in x_pos_0based:
        pi_X[i] = j
        sigma_X[j] = i

    # Step 2: Calculate the writhe of the grid diagram
    writhe = 0
    crossings = []
    
    # Iterate over all possible grid squares to find crossings
    for i in range(n):  # column index i
        for j in range(n):  # row index j
            # A crossing exists if vertical segment in col 'i' and horizontal segment in row 'j' intersect.
            
            # Check for vertical span over row j
            v_span = (pi_O[i] < j < pi_X[i]) or (pi_X[i] < j < pi_O[i])
            
            # Check for horizontal span over column i
            h_span = (sigma_O[j] < i < sigma_X[j]) or (sigma_X[j] < i < sigma_O[j])
            
            if v_span and h_span:
                # A crossing exists at (i, j). Determine its sign.
                # Orientation is X->O (vertical) and O->X (horizontal)
                
                # Vertical orientation: +1 for Up (O is above X), -1 for Down
                v_dir = 1 if pi_O[i] > pi_X[i] else -1
                
                # Horizontal orientation: +1 for Right (X is to the right of O), -1 for Left
                h_dir = 1 if sigma_X[j] > sigma_O[j] else -1
                
                sign = v_dir * h_dir
                writhe += sign
                # Store crossing info with 1-based coordinates for printing
                crossings.append({'pos': (i + 1, j + 1), 'sign': sign})
    
    # Step 3: Calculate the Thurston-Bennequin number using the formula tb = writhe - n
    tb = writhe - n

    # Step 4: Print the detailed calculation
    print("Step 1: Calculating the writhe (w) of the grid diagram.")
    print("The writhe is the sum of the signs of all crossings.")
    
    pos_crossings_count = sum(1 for c in crossings if c['sign'] == 1)
    neg_crossings_count = len(crossings) - pos_crossings_count
    
    print(f"Found {len(crossings)} crossings:")
    for c in sorted(crossings, key=lambda x: x['pos']):
        print(f"  - Crossing at (column={c['pos'][0]}, row={c['pos'][1]}) has sign {c['sign']:+}.")
        
    print(f"\nThe total writhe is the sum of these signs.")
    print(f"w = (Number of +1 crossings) - (Number of -1 crossings)")
    print(f"w = {pos_crossings_count} - {neg_crossings_count} = {writhe}")

    print("\nStep 2: Calculating the Thurston-Bennequin number (tb).")
    print("Using the formula: tb = w - n")
    print(f"Here, n is the grid size, so n = {n}.")
    
    print("\nFinal Calculation:")
    print(f"tb = {writhe} - {n} = {tb}")


solve_thurston_bennequin()
print("<<< -2 >>>")