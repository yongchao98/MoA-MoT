import math

def get_squares_for_segment(p1, p2):
    """
    Numerically determines the set of unique grid squares a line segment passes through.
    The segment is from p1=(x1, y1) to p2=(x2, y2).
    """
    x1, y1 = p1
    x2, y2 = p2
    squares = set()
    
    dx = x2 - x1
    dy = y2 - y1
    
    # The number of steps for the simulation must be large enough to not "jump"
    # over any grid squares. A number of steps slightly larger than the Manhattan
    # distance between the endpoints' integer coordinates is sufficient.
    num_steps = int(abs(dx) + abs(dy)) + 2
    
    for i in range(num_steps + 1):
        t = i / num_steps
        x = x1 + t * dx
        y = y1 + t * dy
        squares.add((math.floor(x), math.floor(y)))
        
    return squares

# --- Define Triangle Geometry in Optimal Position ---

# Use small offsets to shift the triangle off lattice points.
epsilon = 0.001
delta = 0.002

# Define side lengths and projections for the chosen orientation.
side_leg = 18.0
side_hypotenuse = 18.0 * math.sqrt(2)
leg_projection = side_leg / math.sqrt(2)  # This is 9*sqrt(2)

# Define the vertices. The right angle is at C, and the hypotenuse AB is horizontal.
A = (epsilon, delta)
B = (side_hypotenuse + epsilon, delta)
C = (leg_projection + epsilon, leg_projection + delta)

# --- Calculate Squares for Each Side ---

squares_leg_AC = get_squares_for_segment(A, C)
squares_leg_BC = get_squares_for_segment(C, B)
squares_hypotenuse_AB = get_squares_for_segment(A, B)

# --- Final Calculation using Inclusion-Exclusion Principle ---

# Get the size of each set
len_ac = len(squares_leg_AC)
len_bc = len(squares_leg_BC)
len_ab = len(squares_hypotenuse_AB)

# Get the size of overlaps between pairs of sides
intersection_ac_bc = squares_leg_AC.intersection(squares_leg_BC)
intersection_ac_ab = squares_leg_AC.intersection(squares_hypotenuse_AB)
intersection_bc_ab = squares_leg_BC.intersection(squares_hypotenuse_AB)

len_ac_bc = len(intersection_ac_bc)
len_ac_ab = len(intersection_ac_ab)
len_bc_ab = len(intersection_bc_ab)

# Get the size of the overlap of all three sides
intersection_all = intersection_ac_bc.intersection(squares_hypotenuse_AB)
len_all = len(intersection_all)

# The total number of unique squares k is the size of the union of all sets.
k = len(squares_leg_AC.union(squares_leg_BC).union(squares_hypotenuse_AB))
calculated_k = len_ac + len_bc + len_ab - (len_ac_ab + len_bc_ab + len_ac_bc) + len_all

# --- Output the results ---
print("The optimal configuration places the triangle's legs diagonally to the grid.")
print("The total number of squares is calculated using the Principle of Inclusion-Exclusion.")
print("Equation: k = |Squares(AC)| + |Squares(BC)| + |Squares(AB)| - Overlaps")
print(f"k = {len_ac} + {len_bc} + {len_ab} - ({len_ac_ab} + {len_bc_ab} + {len_ac_bc}) + {len_all}")
print(f"k = {calculated_k}")
print(f"\nThe largest number k such that the triangle's perimeter can pass through at least k coordinate grid squares is {k}.")
