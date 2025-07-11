def solve_minimal_critical_points():
    """
    Calculates the minimal number of critical points for a smooth function
    on a 2-torus using Morse theory.
    """
    
    explanation = """
This problem can be solved using a result from Morse theory, which connects the topology of a manifold to the properties of smooth functions defined on it.

A fundamental theorem in this area states that for any smooth function on a compact manifold, the total number of critical points must be greater than or equal to the sum of its Betti numbers. Betti numbers are topological invariants that, roughly speaking, count the number of 'holes' of different dimensions in a space.

For the 2-torus (T^2), the Betti numbers are:
- b_0 = 1 (the torus is one connected piece)
- b_1 = 2 (the torus has two independent loops, one 'around' the donut and one 'through' the hole)
- b_2 = 1 (the torus encloses one 2D 'void')

The minimal number of critical points is therefore the sum of these Betti numbers. This bound is known to be sharp, as there exist functions on the torus that have exactly this many critical points.
"""
    print(explanation)

    # Betti numbers for the 2-torus
    b_0 = 1
    b_1 = 2
    b_2 = 1
    
    # Calculate the minimal number of critical points by summing the Betti numbers
    min_points = b_0 + b_1 + b_2
    
    # Output the final calculation as an equation
    print("The calculation is as follows:")
    print(f"Minimal Number of Critical Points = b_0 + b_1 + b_2 = {b_0} + {b_1} + {b_2} = {min_points}")

solve_minimal_critical_points()
