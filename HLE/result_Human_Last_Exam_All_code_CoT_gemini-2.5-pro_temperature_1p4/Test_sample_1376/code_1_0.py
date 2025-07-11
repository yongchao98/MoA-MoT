def solve_tomb_engraving():
    """
    Calculates the optimal number of squares and circles to maximize
    the number of engraved characters on a meteorite sheet.
    """
    # Define constants based on the problem description
    sheet_w = 140  # cm
    sheet_h = 110  # cm
    square_size = 10  # cm
    circle_diameter = 40  # cm (from 20cm radius)
    chars_per_square = 4
    chars_per_circle = 999

    # The greedy strategy is to maximize circles first, as they are more valuable per area.
    # We need to check both orientations of the sheet to find the global maximum.
    
    best_config = {'N': 0, 'M': 0, 'K': 0}

    # Iterate through two possible orientations of the sheet
    orientations = [(sheet_w, sheet_h), (sheet_h, sheet_w)]
    
    for i, (w, h) in enumerate(orientations):
        # Calculate how many 40x40 blocks for circles can fit in a grid
        m_fit_w = w // circle_diameter
        m_fit_h = h // circle_diameter
        M = m_fit_w * m_fit_h

        # Calculate the area consumed by the grid of circles
        used_w = m_fit_w * circle_diameter
        used_h = m_fit_h * circle_diameter

        # The remaining area is an L-shape. We split it into two rectangles to
        # calculate the number of squares that can fit.
        # This split gives two rectangles: (w - used_w) x h and used_w x (h - used_h)
        rem_rect1_w = w - used_w
        rem_rect1_h = h
        
        rem_rect2_w = used_w
        rem_rect2_h = h - used_h

        # Calculate how many 10x10 squares fit in each remaining rectangle
        n_in_rect1 = (rem_rect1_w // square_size) * (rem_rect1_h // square_size)
        n_in_rect2 = (rem_rect2_w // square_size) * (rem_rect2_h // square_size)
        N = n_in_rect1 + n_in_rect2

        # Calculate total characters for this configuration
        K = (N * chars_per_square) + (M * chars_per_circle)
        
        # If this configuration is the best so far, save it
        if K > best_config['K']:
            best_config['N'] = N
            best_config['M'] = M
            best_config['K'] = K

    # Unpack the best results
    N = best_config['N']
    M = best_config['M']
    K = best_config['K']

    # Print the final results
    print(f"To maximize the number of engraved characters, the workers should produce:")
    print(f"- {N} squares (N)")
    print(f"- {M} circles (M)")
    print("\nThis is based on the following calculation:")
    print(f"({N} squares * {chars_per_square} chars) + ({M} circles * {chars_per_circle} chars) = {K} total characters")
    
    final_answer = f"{N}:{M}:{K}"
    print(f"\nThe final answer in the format N:M:K is: {final_answer}")
    
    # Do not remove the following line
    print(f"\n<<<{final_answer}>>>")

solve_tomb_engraving()