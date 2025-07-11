import math

def get_volume_set_description(d):
    """
    Determines the set of possible values for the volume of the polytope P
    constructed from a d-dimensional simplex T with volume V.

    Args:
        d (int): The dimension of the simplex.
    """
    if not isinstance(d, int) or d < 1:
        print("Error: Dimension d must be a positive integer.")
        return

    # A d-dimensional simplex has d+1 vertices and C(d+1, 2) edges.
    # The problem statement's formula n = C(d, 2) seems to be a typo.
    num_edges = math.comb(d + 1, 2)
    
    # The polytope P is bounded by 2n hyperplanes.
    num_hyperplanes = 2 * num_edges

    fact_d = math.factorial(d)
    V_str = "V"

    print(f"For a simplex in dimension d = {d}:")
    print(f"The polytope P is bounded by {num_hyperplanes} hyperplanes.")

    if d == 1:
        print("The set of possible values for the volume of P is a single point.")
        print("The final equation for the volume of P is:")
        print(f"Volume(P) = {fact_d} * {V_str}")
    elif d == 2:
        print("The set of possible values for the volume of P is a single point.")
        print("This is because all triangles (2D simplices) are orthocentric.")
        print("The final equation for the volume of P is:")
        print(f"Volume(P) = {fact_d} * {V_str}")
    elif d == 3:
        lower_bound = fact_d
        upper_bound = 8
        print("The set of possible values for the volume of P is a closed interval.")
        print(f"The minimum value, {lower_bound}*V, is achieved for orthocentric tetrahedra.")
        print(f"The maximum value is {upper_bound}*V.")
        print("The final relation for the volume of P is:")
        print(f"{lower_bound} * {V_str} <= Volume(P) <= {upper_bound} * {V_str}")
    else:  # d >= 4
        lower_bound = fact_d
        print("The set of possible values for the volume of P is an unbounded interval.")
        print(f"The minimum value, {lower_bound}*V, is achieved for orthocentric simplices.")
        print("The volume is not bounded from above.")
        print("The final relation for the volume of P is:")
        print(f"Volume(P) >= {lower_bound} * {V_str}")

if __name__ == '__main__':
    # Demonstrate for several dimensions
    get_volume_set_description(1)
    print("-" * 20)
    get_volume_set_description(2)
    print("-" * 20)
    get_volume_set_description(3)
    print("-" * 20)
    get_volume_set_description(4)
