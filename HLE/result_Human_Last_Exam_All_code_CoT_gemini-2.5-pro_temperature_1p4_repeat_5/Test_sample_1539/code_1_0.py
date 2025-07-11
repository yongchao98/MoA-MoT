import math

def solve_book_engraving():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters, and prints the result.
    """

    # Step 1: Define constants from the problem description
    
    # Material dimensions
    rect_width_cm = 140
    rect_height_cm = 110

    # Circle ("Book of Heaven") parameters
    circle_radius_cm = 20
    circle_diameter_cm = 2 * circle_radius_cm
    circle_symbols = 9999

    # Square ("Book of Earth") parameters
    square_side_cm = 10
    square_symbols = 360

    # Encoding parameters
    num_chars_to_encode = 1000
    yinyang_wuxing_alphabets = 10  # 2 (yinyang) * 5 (wuxing)
    bagua_alphabets = 8

    # Step 2: Calculate the number of characters per plate
    
    # For circular plates
    symbols_per_char_circle = math.ceil(math.log(num_chars_to_encode, yinyang_wuxing_alphabets))
    chars_per_circle = math.floor(circle_symbols / symbols_per_char_circle)

    # For square plates
    symbols_per_char_square = math.ceil(math.log(num_chars_to_encode, bagua_alphabets))
    chars_per_square = math.floor(square_symbols / symbols_per_char_square)

    # Step 3: Solve the packing problem
    # Due to much higher value density (chars/cm^2), we prioritize packing circles.
    # We use a hexagonal packing strategy for maximum density.

    # With width=140, height=110:
    # Height for 3 staggered rows = D + 2*(R*sqrt(3)) = 40 + 2*(20*sqrt(3)) approx 109.28cm < 110cm. So 3 rows fit.
    # Row 1 (non-staggered): floor(140/40) = 3 circles
    # Row 2 (staggered): floor((140-R)/D) = floor((140-20)/40) = 3 circles
    # Row 3 (non-staggered): floor(140/40) = 3 circles
    # This gives a total of 3+3+3 = 9 circles.
    num_circles = 9

    # The leftover space from this highly efficient packing is in small, curved
    # cusps between circles, which are too small to yield a 10x10cm square.
    num_squares = 0

    # Step 4: Calculate the total maximum number of characters
    max_chars = (num_circles * chars_per_circle) + (num_squares * chars_per_square)

    # Step 5: Print the final result in the required format
    print("Final Plan Calculation:")
    print(f"Number of circular plates (N): {num_circles}")
    print(f"Number of square plates (M): {num_squares}")
    print("\nEquation for total characters (K):")
    # The user requested to output each number in the final equation
    print(f"K = {num_circles} * {chars_per_circle} + {num_squares} * {chars_per_square} = {max_chars}")
    
    print("\nFinal Answer Format (N:M:K):")
    final_answer = f"{num_circles}:{num_squares}:{max_chars}"
    print(final_answer)
    return final_answer

# Execute the function
final_answer_string = solve_book_engraving()

# The final answer in the required format for submission
# print(f"<<<{final_answer_string}>>>")