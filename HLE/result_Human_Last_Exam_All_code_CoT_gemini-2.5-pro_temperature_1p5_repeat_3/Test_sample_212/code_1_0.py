import sys

# We are asked if two actions of the fundamental group on the fiber of a universal cover are the same for the torus T^2.
# This script will simulate both actions on an arbitrary point in the fiber to demonstrate the answer.

# In our case, X = T^2 (the torus).
# The universal cover X_tilde is R^2 (the plane).
# The fundamental group pi_1(T^2) is Z x Z. We denote its generators by 'alpha' and 'beta'.
# The fiber p_inv(x0) over the basepoint is the integer lattice Z^2.

class Point:
    """A point in the fiber Z^2, which is a subset of R^2."""
    def __init__(self, m, n):
        self.m = int(m)
        self.n = int(n)
        
    def __eq__(self, other):
        return self.m == other.m and self.n == other.n
        
    def __repr__(self):
        return f"({self.m}, {self.n})"

def action_by_holonomy(loop_generator, point):
    """
    Computes the action by holonomy (path lifting).
    The action of a loop on a point is the endpoint of the lift of that loop starting from the point.
    - The loop 'alpha' in T^2 can be represented by the path gamma(t) = (t, 0) for t in [0,1].
    - Lifting this path from a starting point (m, n) in R^2 gives the path gamma_tilde(t) = (m + t, n).
    - The endpoint at t=1 is (m + 1, n).
    - Similarly, the lift of 'beta' from (m,n) ends at (m, n+1).
    """
    if loop_generator == 'alpha':
        return Point(point.m + 1, point.n)
    elif loop_generator == 'beta':
        return Point(point.m, point.n + 1)
    else:
        raise ValueError("Unknown loop generator.")

def action_by_deck_transformation(loop_generator, point):
    """
    Computes the action by deck transformations.
    First, find the deck transformation corresponding to the loop. This is done by lifting the
    loop from the fiber's basepoint (0,0). The resulting endpoint defines the translation vector
    of the deck transformation.
    """
    base_point_in_fiber = Point(0, 0)
    
    # The endpoint of the lift of the loop from the basepoint gives the translation vector.
    # This is equivalent to applying the holonomy action to the base point.
    endpoint = action_by_holonomy(loop_generator, base_point_in_fiber)
    translation_vector = (endpoint.m, endpoint.n)
    
    # The action is to apply this deck transformation (translation) to the point.
    return Point(point.m + translation_vector[0], point.n + translation_vector[1])

def solve():
    """
    Performs the comparison and prints the final answer.
    """
    # Let's choose an arbitrary point in the fiber Z^2 to test the actions.
    # We take command line arguments if provided, otherwise default to (5, 3).
    m = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    test_point = Point(m, n)
    
    print(f"Let's test the actions on an arbitrary point in the fiber: {test_point}")
    print("-" * 50)
    
    # --- Test for generator 'alpha' ---
    loop_alpha = 'alpha'
    print(f"Action of generator '{loop_alpha}':")
    
    holonomy_result_alpha = action_by_holonomy(loop_alpha, test_point)
    print(f"1. Action by holonomy:")
    print(f"   The lift of loop '{loop_alpha}' starting at {test_point} ends at {holonomy_result_alpha}.")
    
    deck_result_alpha = action_by_deck_transformation(loop_alpha, test_point)
    print(f"2. Action by deck transformation:")
    print(f"   The deck transformation for '{loop_alpha}' is translation by (1, 0).")
    print(f"   Applying this to {test_point} gives {deck_result_alpha}.")

    print(f"   Results for '{loop_alpha}' are the same: {holonomy_result_alpha == deck_result_alpha}")
    print("-" * 50)

    # --- Test for generator 'beta' ---
    loop_beta = 'beta'
    print(f"Action of generator '{loop_beta}':")
    
    holonomy_result_beta = action_by_holonomy(loop_beta, test_point)
    print(f"1. Action by holonomy:")
    print(f"   The lift of loop '{loop_beta}' starting at {test_point} ends at {holonomy_result_beta}.")
    
    deck_result_beta = action_by_deck_transformation(loop_beta, test_point)
    print(f"2. Action by deck transformation:")
    print(f"   The deck transformation for '{loop_beta}' is translation by (0, 1).")
    print(f"   Applying this to {test_point} gives {deck_result_beta}.")
    
    print(f"   Results for '{loop_beta}' are the same: {holonomy_result_beta == deck_result_beta}")
    print("-" * 50)
    
    print("Conclusion: The two actions produce identical results for the generators of the fundamental group.")
    print("Since the actions agree on the generators, they agree for all elements of the group.")

if __name__ == "__main__":
    solve()
    print("\nAre these two actions the same when X = T^2?")
    print("<<<Yes>>>")
