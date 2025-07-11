def solve_critical_points_torus():
    """
    Calculates the minimal number of critical points for a smooth function
    f: T^2 -> R using Morse theory.
    """

    # Step 1: Define the Betti numbers for the 2-torus (T^2).
    # The Betti numbers are topological invariants of a space.
    # b_0: number of connected components. The torus is one piece.
    b0 = 1
    # b_1: number of 1-dimensional "holes". The torus has two fundamental loops.
    b1 = 2
    # b_2: number of 2-dimensional "voids". The torus encloses one void.
    b2 = 1

    print("Step 1: The Betti numbers for the 2-torus (T^2) are:")
    print(f"b_0 (connected components) = {b0}")
    print(f"b_1 (1D holes/loops)      = {b1}")
    print(f"b_2 (2D voids)            = {b2}")
    print("-" * 50)

    # Step 2: Apply the Morse inequalities.
    # For a Morse function, the number of critical points of index k (c_k) is
    # greater than or equal to the k-th Betti number (b_k).
    # c_0 are local minima, c_1 are saddle points, c_2 are local maxima.
    print("Step 2: According to the Morse inequalities (c_k >= b_k), the minimal number of critical points of each type is:")
    print(f"Minimal number of minima (c_0)  >= b_0  => c_0_min = {b0}")
    print(f"Minimal number of saddles (c_1) >= b_1  => c_1_min = {b1}")
    print(f"Minimal number of maxima (c_2)  >= b_2  => c_2_min = {b2}")
    print("-" * 50)

    # Step 3: Calculate the minimal total number of critical points.
    # A theorem states that the total number of critical points for any smooth function
    # is at least the sum of the Betti numbers.
    min_total_critical_points = b0 + b1 + b2

    print("Step 3: Calculate the minimal total number of critical points.")
    print("This is the sum of the minimal numbers for each index, which corresponds to the sum of the Betti numbers.")
    print("The final equation is:")
    print(f"Minimal Total = {b0} + {b1} + {b2} = {min_total_critical_points}")
    print("-" * 50)
    
    # Final check with Euler Characteristic: sum((-1)^k * c_k) = sum((-1)^k * b_k)
    # Our choice of c_k = b_k satisfies this trivially: b0 - b1 + b2 = b0 - b1 + b2
    euler_char = b0 - b1 + b2
    print(f"This result is consistent with the Euler characteristic of the torus, chi(T^2) = {b0} - {b1} + {b2} = {euler_char}.")
    
    print("\nSince a function with exactly this number of critical points exists, this is the minimal number.")
    print(f"\nFinal Answer: The minimal number of critical points is {min_total_critical_points}.")

if __name__ == '__main__':
    solve_critical_points_torus()