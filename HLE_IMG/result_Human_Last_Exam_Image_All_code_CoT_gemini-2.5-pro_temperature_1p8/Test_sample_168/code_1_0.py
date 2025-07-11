import itertools

def main():
    """
    Calculates the number of ways a light ray can travel from M to N with
    2 reflections on G1 and 1 on G2, G3, and G4.
    """
    # Step 1: Find all physically possible reflection sequences.
    vertical_mirrors = ['G1', 'G1', 'G3']
    horizontal_mirrors = ['G2', 'G4']

    # Get unique permutations for V and H mirror sets
    v_perms = set(itertools.permutations(vertical_mirrors))
    h_perms = set(itertools.permutations(horizontal_mirrors))
    
    # Total sequences follow the pattern V-H-V-H-V
    sequences = []
    for v_perm in v_perms:
        for h_perm in h_perms:
            # Weave them together: V1, H1, V2, H2, V3
            seq = [v_perm[0], h_perm[0], v_perm[1], h_perm[1], v_perm[2]]
            sequences.append(seq)
    
    num_sequences = len(sequences)
    print(f"The number of ways to arrange the 3 vertical mirrors (G1, G1, G3) is {len(v_perms)}.")
    print(f"The number of ways to arrange the 2 horizontal mirrors (G2, G4) is {len(h_perms)}.")
    print(f"Total number of possible reflection sequences: {len(v_perms)} * {len(h_perms)} = {num_sequences}\n")
    
    # Step 2 & 3: Model reflections and define transformations.
    # We represent a virtual image's coordinates (x', y') as coefficients:
    # x' = ax*x + aW*W, y' = ay*y + aH*H
    # So a point is a tuple of coefficients (ax, aW, ay, aH).
    # N is initially (x, y), so its coefficient representation is (1, 0, 1, 0).
    def t1(p): return (-p[0], -p[1], p[2], p[3])
    def t2(p): return (p[0], p[1], -p[2], -p[3])
    def t3(p): return (-p[0], 2 - p[1], p[2], p[3])
    def t4(p): return (p[0], p[1], -p[2], 2 - p[3])

    mirror_ops = {'G1': t1, 'G2': t2, 'G3': t3, 'G4': t4}
    
    # Step 4: Calculate unique final images.
    final_images = set()

    # The light path M -> R1 -> R2 -> ... -> N corresponds to a straight line
    # from M to the virtual image N' = T_R1(T_R2(...(T_Rk(N))...)).
    # To compute this, we apply the transformations for Rk, Rk-1, ..., R1 in order.
    # This corresponds to iterating through the sequence in reverse.
    print("Calculating the virtual image N' for each sequence:")
    for i, seq in enumerate(sequences):
        current_image = (1, 0, 1, 0)  # Start with N = (x, y)
        for mirror_name in reversed(seq):
            op = mirror_ops[mirror_name]
            current_image = op(current_image)
        final_images.add(current_image)
        print(f"Sequence {i+1}: {seq} -> Image Coeffs: {current_image}")
        
    num_ways = len(final_images)

    print("\n--- Final Results ---")
    print(f"Out of {num_sequences} possible sequences, there are {num_ways} unique virtual images.")
    print("Each unique virtual image corresponds to a distinct light ray path.")
    
    print("\nThe coordinates of the unique virtual images (x', y') are:")
    for i, image_coeffs in enumerate(sorted(list(final_images))):
        ax, aW, ay, aH = image_coeffs
        print(f"Path {i+1}: x' = {ax}*x + {aW}*W, y' = {ay}*y + {aH}*H")

    print(f"\nThus, the total number of ways to draw the light ray path is {num_ways}.")
    print(f"\n<<<4>>>")

if __name__ == '__main__':
    main()