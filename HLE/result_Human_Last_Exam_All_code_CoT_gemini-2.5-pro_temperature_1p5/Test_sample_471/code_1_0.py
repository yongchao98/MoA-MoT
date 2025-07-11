def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """
    
    # The Betti numbers for the 2-torus (T^2) are as follows:
    # b_0: number of connected components
    b0 = 1
    # b_1: number of 1-dimensional "circular" holes
    b1 = 2
    # b_2: number of 2-dimensional "cavities"
    b2 = 1
    
    # According to Morse theory, the minimal number of critical points of a
    # smooth function on a compact manifold is the sum of its Betti numbers.
    # A simple example function on the torus that achieves this minimum has:
    # 1 minimum (index 0, corresponding to b0)
    # 2 saddle points (index 1, corresponding to b1)
    # 1 maximum (index 2, corresponding to b2)
    min_critical_points = b0 + b1 + b2
    
    # Print the equation that demonstrates this calculation.
    print(f"{b0} + {b1} + {b2} = {min_critical_points}")

solve_minimal_critical_points()