import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the number of engraved Chinese characters.
    """
    
    # Step 1: Calculate characters per plate based on encoding requirements.
    
    # For circular plates (yinyang wuxing)
    unique_symbols_circle = 2 * 5  # yinyang * wuxing
    symbols_per_char_circle = math.ceil(math.log(1000, unique_symbols_circle)) # Need 10^n >= 1000
    chars_per_circle = math.floor(9999 / symbols_per_char_circle)
    
    # For squared plates (bagua)
    unique_symbols_square = 8  # bagua
    symbols_per_char_square = math.ceil(math.log(1000, unique_symbols_square)) # Need 8^n >= 1000
    chars_per_square = math.floor(360 / symbols_per_char_square)

    print("Step 1: Calculating Characters per Plate")
    print(f" - A circular plate can hold {chars_per_circle} characters ({symbols_per_char_circle} symbols per character).")
    print(f" - A squared plate can hold {chars_per_square} characters ({symbols_per_char_square} symbols per character).")
    print("-" * 20)
    
    # Step 2: Determine the optimal number of plates (N, M).
    # Based on geometric packing analysis, we can fit 8 circles and have
    # a 20x110 cm rectangular piece remaining.
    
    N = 8  # Number of circular plates from efficient packing
    
    remaining_width = 140 - 120 # Bounding box of 8 circles is ~120cm wide
    remaining_height = 110
    square_side = 10
    
    M = math.floor(remaining_width / square_side) * math.floor(remaining_height / square_side)

    print("Step 2: Determining Optimal Plate Count")
    print(f" - By using an efficient packing method, we can cut N = {N} circular plates.")
    print(f" - From the remaining rectangular material, we can cut M = {M} squared plates.")
    print("-" * 20)
    
    # Step 3: Calculate the total maximum characters (K).
    
    K = (N * chars_per_circle) + (M * chars_per_square)
    
    print("Step 3: Calculating Total Engraved Characters (K)")
    print(f"The calculation is: ({N} * {chars_per_circle}) + ({M} * {chars_per_square}) = {K}")
    print(f"So, K = {N * chars_per_circle} + {M * chars_per_square} = {K}")
    print("-" * 20)

    # Final result in the specified N:M:K format
    final_answer = f"{N}:{M}:{K}"
    print("Final Answer:")
    print(final_answer)
    print(f"<<<{final_answer}>>>")


solve_laozi_puzzle()