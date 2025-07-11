import math

def solve_emperor_tomb_problem():
    """
    Solves the optimization problem for Emperor Qin's artifacts.
    This function calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Step 1: Define the parameters from the problem description.
    sheet_w = 140
    sheet_h = 110
    square_dim = 10
    circle_radius = 20
    circle_dim = circle_radius * 2  # Bounding box for the circle is its diameter

    chars_per_square = 4  # "Qin Shi Huang Di"
    symbols_per_circle = 999
    
    # --- Calculation for Characters per Circle ---
    # The bio has a set of 7x7x9 unique characters.
    unique_bio_chars = 7 * 7 * 9
    # The bio is encoded using 8 unique Bagua symbols.
    # We need to find how many symbols make one character.
    # This is equivalent to finding 'x' in 8^x >= unique_bio_chars
    symbols_per_char = math.ceil(math.log(unique_bio_chars, 8))
    # Calculate how many characters fit on one circle.
    chars_per_circle = math.floor(symbols_per_circle / symbols_per_char)
    
    print("--- Step 1: Problem Setup ---")
    print(f"The meteorite material is a {sheet_w}x{sheet_h} cm rectangle.")
    print(f"Squares are {square_dim}x{square_dim} cm and hold {chars_per_square} characters each.")
    print(f"Circles need a {circle_dim}x{circle_dim} cm area and are engraved with symbols.")
    print(f"The bio has {unique_bio_chars} unique characters, requiring {symbols_per_char} Bagua symbols per character.")
    print(f"A circle holds {symbols_per_circle} symbols, so it can contain {symbols_per_circle} / {symbols_per_char} = {chars_per_circle} characters.\n")
    
    # Step 2: Formulate the optimization problem.
    # K = 4*N + 333*M
    # The constraint is packing. Since all dimensions are multiples of 10,
    # we can discretize the problem into a grid of 10x10 cells.
    total_cells = (sheet_w / square_dim) * (sheet_h / square_dim)
    cells_per_circle = (circle_dim / square_dim) * (circle_dim / square_dim)
    
    print("--- Step 2: Optimization Approach ---")
    print("The goal is to maximize K = 4*N + 333*M.")
    print(f"The {sheet_w}x{sheet_h} rectangle can be seen as a grid of {int(total_cells)} cells of size 10x10.")
    print(f"A square occupies 1 cell. A circle occupies a {int(circle_dim/10)}x{int(circle_dim/10)} area, which is {int(cells_per_circle)} cells.")
    print("The packing constraint simplifies to: N + 16*M = 154.")
    print("Substituting N into the goal equation: K = 4*(154 - 16*M) + 333*M = 616 + 269*M.")
    print("To maximize K, we must maximize M.\n")
    
    # Step 3: Find the maximum possible value for M.
    max_m_orientation1 = math.floor(sheet_w / circle_dim) * math.floor(sheet_h / circle_dim)
    # Check the other orientation of the rectangle (110x140)
    max_m_orientation2 = math.floor(sheet_h / circle_dim) * math.floor(sheet_w / circle_dim)
    max_m = max(max_m_orientation1, max_m_orientation2)
    
    print("--- Step 3: Find the Maximum Number of Circles (M) ---")
    print(f"Maximum circles that fit along the 140cm side: floor(140 / 40) = {math.floor(sheet_w / circle_dim)}")
    print(f"Maximum circles that fit along the 110cm side: floor(110 / 40) = {math.floor(sheet_h / circle_dim)}")
    print(f"Total possible circles M = {math.floor(sheet_w / circle_dim)} * {math.floor(sheet_h / circle_dim)} = {max_m}.\n")

    # Step 4: Calculate the final numbers for N, M, and K.
    # The optimal solution is at M = max_m.
    optimal_m = max_m
    optimal_n = total_cells - (optimal_m * cells_per_circle)
    optimal_k = (optimal_n * chars_per_square) + (optimal_m * chars_per_circle)
    
    print("--- Step 4: Calculate the Optimal Solution ---")
    print(f"Using the maximum M = {optimal_m}:")
    print(f"Number of squares N = {int(total_cells)} - {int(cells_per_circle)} * {optimal_m} = {int(optimal_n)}")
    print("The maximum number of characters K is calculated as:")
    print(f"K = (N * chars_per_square) + (M * chars_per_circle)")
    print(f"K = ({int(optimal_n)} * {chars_per_square}) + ({optimal_m} * {chars_per_circle}) = {int(optimal_n * chars_per_square)} + {int(optimal_m * chars_per_circle)} = {int(optimal_k)}\n")
    
    # Step 5: Final Answer
    print("--- Final Answer ---")
    print(f"The optimal combination is N={int(optimal_n)} squares and M={optimal_m} circles.")
    print(f"This yields a maximum of K={int(optimal_k)} engraved characters.")
    print("The final answer in the format N:M:K is:")
    print(f"{int(optimal_n)}:{optimal_m}:{int(optimal_k)}")

# Run the solver
solve_emperor_tomb_problem()