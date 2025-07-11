import sympy

def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function on a 2-torus
    using Morse theory.
    """
    
    print("Step 1: Define the Betti numbers for the 2-torus (T^2).")
    # b_k is the k-th Betti number.
    # For a 2-torus T^2:
    # b_0 = number of connected components = 1
    # b_1 = number of 1D holes (loops) = 2
    # b_2 = number of 2D holes (voids) = 1
    betti_numbers = {
        0: 1,
        1: 2,
        2: 1
    }
    print(f"b_0 = {betti_numbers[0]} (minima)")
    print(f"b_1 = {betti_numbers[1]} (saddles)")
    print(f"b_2 = {betti_numbers[2]} (maxima)")
    print("-" * 30)

    print("Step 2: Apply the Morse Inequalities.")
    print("The number of critical points of index k (c_k) must be at least the k-th Betti number (b_k).")
    print("c_k >= b_k")
    print("-" * 30)

    print("Step 3: Calculate the minimal number of critical points for each index.")
    # c_k_min is the minimal number of critical points of index k.
    c_min = {}
    for k, b_k in betti_numbers.items():
        c_min[k] = b_k # The minimal number is the Betti number itself.
        print(f"Minimal c_{k} = b_{k} = {b_k}")
    
    print("-" * 30)
    
    print("Step 4: Calculate the total minimal number of critical points.")
    total_min_critical_points = sum(c_min.values())
    
    # Create the equation string
    equation_parts = [f"{c}" for c in c_min.values()]
    equation_str = " + ".join(equation_parts)
    
    print(f"The minimal total number of critical points is the sum of the minimal c_k values.")
    print(f"Minimal Total = c_0 + c_1 + c_2 = {equation_str} = {total_min_critical_points}")
    print("-" * 30)
    
    # Bonus: Check the Euler characteristic
    # The alternating sum of critical points must equal the Euler characteristic.
    # Euler characteristic chi(T^2) = b_0 - b_1 + b_2 = 1 - 2 + 1 = 0
    # For our minimal numbers: c_0 - c_1 + c_2 = 1 - 2 + 1 = 0. It holds.
    
    print(f"The minimal number of critical points for a smooth function f: T^2 -> R is {total_min_critical_points}.")
    

solve_minimal_critical_points()