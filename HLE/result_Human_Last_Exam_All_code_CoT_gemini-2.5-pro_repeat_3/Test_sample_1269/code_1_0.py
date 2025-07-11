import sys

def solve():
    """
    Solves the Hfr mapping problem by finding the closest landmark gene to the target gene.
    """
    # Step 1: Define the E. coli genetic map (positions in minutes)
    gene_map = {
        'thr': 0,
        'azis': 2,
        'pro': 6,
        'lac': 8,
        'ton': 27,
        'str': 73
    }
    
    # Step 2: Identify the target gene and the landmark genes from the options
    target_gene = 'azis'
    landmarks = ['ton', 'lac', 'pro', 'thr', 'str']
    
    print("Plan:")
    print("1. The problem states the 'azis' gene is transferred very early.")
    print("2. This means the origin of transfer (O) must be located right next to 'azis' on the chromosome.")
    print(f"3. We will find which landmark gene from the options is closest to 'azis' (at {gene_map[target_gene]} min) on the 100-minute circular map.")
    print("-" * 20)
    
    target_pos = gene_map[target_gene]
    
    min_distance = float('inf')
    closest_landmark = None
    
    print("Calculating distances to 'azis':")
    
    # Step 3: Calculate the distance from 'azis' to each landmark
    for landmark in landmarks:
        landmark_pos = gene_map[landmark]
        
        # Distance on a circle is the minimum of the direct path and the path wrapping around
        distance = abs(target_pos - landmark_pos)
        circular_distance = min(distance, 100 - distance)
        
        print(f" - Distance from '{target_gene}' ({target_pos} min) to '{landmark}' ({landmark_pos} min) is {circular_distance} minutes.")
        
        if circular_distance < min_distance:
            min_distance = circular_distance
            closest_landmark = landmark
            
    print("-" * 20)
    print(f"Conclusion: The closest landmark to '{target_gene}' is '{closest_landmark}'.")
    print("Therefore, the origin of transfer must be near 'thr'.")
    print("This corresponds to the option describing an origin near 'thr'.")
    print("\nFinal Answer Choice Text:")
    print("D. Counterclockwise direction, origin near thr")

solve()
<<<D>>>