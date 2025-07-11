import numpy as np

def main():
    """
    This script demonstrates that a line segment of length 2 can exist on the
    unit sphere of a normed vector space, specifically using the L1-norm in R^2.
    This confirms that the maximum possible length is 2.
    """

    # Define a function for the L1 norm (Manhattan distance)
    def l1_norm(v):
        return np.sum(np.abs(v))

    # Define two points, a and b, in R^2
    a = np.array([1, 0])
    b = np.array([0, 1])

    print(f"Consider the space R^2 with the L1-norm: ||(x,y)||_1 = |x| + |y|.")
    print(f"Let point a = {a} and point b = {b}.")
    print("-" * 40)

    # Verify that a and b lie on the L1 unit sphere
    norm_a = l1_norm(a)
    norm_b = l1_norm(b)
    print(f"The norm of a is ||{a}||_1 = |{a[0]}| + |{a[1]}| = {norm_a}.")
    print(f"The norm of b is ||{b}||_1 = |{b[0]}| + |{b[1]}| = {norm_b}.")
    print("Since their norms are 1, both points are on the unit sphere.")
    print("-" * 40)

    # The line segment from a to b is composed of points x(t) = (1-t)a + tb.
    # For t in [0,1], x(t) = (1-t, t).
    # The norm is ||x(t)||_1 = |1-t| + |t| = (1-t) + t = 1.
    # This shows the entire segment is on the unit sphere.

    # Calculate the length of the segment, which is the norm of the difference vector.
    diff_vector = a - b
    length = l1_norm(diff_vector)

    print("The length of the line segment is ||a - b||_1.")
    print(f"The difference vector is a - b = {a} - {b} = {diff_vector}.")
    print("The length calculation is:")
    # Output each number in the final equation
    print(f"||{diff_vector}||_1 = |{diff_vector[0]}| + |{diff_vector[1]}| = {np.abs(diff_vector[0])} + {np.abs(diff_vector[1])} = {length}")
    print("-" * 40)
    
    print(f"The largest possible length is {length}.")

if __name__ == "__main__":
    main()