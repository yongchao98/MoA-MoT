import math

def solve_emperor_tomb_problem():
    """
    Solves the optimization problem for engraving Qin Shi Huang's artifacts.
    """
    # 1. Define constants from the problem statement
    sheet_w, sheet_h = 140, 110
    square_dim = 10
    circle_radius = 20
    circle_dim = circle_radius * 2

    # 2. Calculate characters per artifact
    chars_per_square = 4
    
    unique_bio_chars = 7 * 7 * 9
    bagua_symbols_count = 8
    symbols_per_circle = 999
    
    # Each character needs to be encoded. We find the number of symbols needed per character.
    # symbols_per_char = ceil(log_base_8(441))
    symbols_per_char = math.ceil(math.log(unique_bio_chars, bagua_symbols_count))
    
    # Calculate how many characters can be fit onto a single circle
    chars_per_circle = symbols_per_circle // symbols_per_char
    
    print("Step 1: Calculating the number of characters per artifact.")
    print(f" - Each square for the name has {chars_per_square} characters.")
    print(f" - The bio uses {unique_bio_chars} unique characters.")
    print(f" - Encoding with {bagua_symbols_count} Bagua symbols requires {symbols_per_char} symbols per character.")
    print(f" - A circle holds {symbols_per_circle} symbols, so it can hold {symbols_per_circle} / {symbols_per_char} = {chars_per_circle} characters.\n")
    
    # 3. Analyze value to prioritize packing circles over squares
    area_for_circle = circle_dim * circle_dim
    squares_in_circle_area = (circle_dim // square_dim) ** 2
    value_of_circle = chars_per_circle
    value_of_squares_in_same_area = squares_in_circle_area * chars_per_square
    
    print("Step 2: Determining the optimization strategy.")
    print(f"A {circle_dim}x{circle_dim}cm area can be used for:")
    print(f" - 1 circle, yielding {value_of_circle} characters.")
    print(f" - {squares_in_circle_area} squares, yielding {squares_in_circle_area} * {chars_per_square} = {value_of_squares_in_same_area} characters.")
    print(f"Since {value_of_circle} > {value_of_squares_in_same_area}, we will maximize the number of circles first.\n")

    # 4. Calculate maximum number of circles (M)
    # Check both orientations of the sheet, although in this case they yield the same M.
    m_option1 = (sheet_w // circle_dim) * (sheet_h // circle_dim)
    m_option2 = (sheet_h // circle_dim) * (sheet_w // circle_dim)
    M = max(m_option1, m_option2)
    
    # 5. Calculate remaining area and the number of squares (N)
    # The configuration for M=6 is a 3x2 grid of 40x40cm blocks
    used_width = (sheet_w // circle_dim) * circle_dim
    used_height = (sheet_h // circle_dim) * circle_dim
    
    total_area = sheet_w * sheet_h
    used_area = used_width * used_height
    remaining_area = total_area - used_area
    
    square_area = square_dim * square_dim
    N = remaining_area // square_area

    print("Step 3: Calculating the optimal number of circles (M) and squares (N).")
    print(f"The 140x110cm sheet can fit a maximum of {M} circles (in a 3x2 grid).")
    print(f"This uses an area of {used_width}x{used_height}cm.")
    print(f"Remaining area is {total_area} - {used_area} = {remaining_area} cm^2.")
    print(f"This remaining area can be cut into {N} squares of {square_dim}x{square_dim}cm.\n")

    # 6. Calculate the maximum total characters (K)
    K = (chars_per_square * N) + (chars_per_circle * M)
    
    print("Step 4: Calculating the maximum total characters (K).")
    print(f"The final optimal numbers are N = {N} squares and M = {M} circles.")
    print("The final equation for the total number of characters K is:")
    print(f"K = (chars on squares) + (chars on circles)")
    print(f"K = ({chars_per_square} * {N}) + ({chars_per_circle} * {M})")
    print(f"K = {chars_per_square * N} + {chars_per_circle * M}")
    print(f"K = {K}\n")

    # 7. Print the final answer in the required format
    final_answer = f"{N}:{M}:{K}"
    print(f"Final Answer (N:M:K) = {final_answer}")

solve_emperor_tomb_problem()