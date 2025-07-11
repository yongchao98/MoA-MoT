import math

def main():
    """
    This script demonstrates that the maximum length of a line segment
    on a unit sphere can be 2, using the L1-norm as an example.
    """

    # Define the L1-norm (Manhattan norm)
    def l1_norm(vector):
        return sum(abs(c) for c in vector)

    # 1. Define two points, u and v, on the L1 unit sphere in R^2.
    u = (1, 0)
    v = (0, 1)

    print(f"Let's consider the L1-norm: ||(x,y)|| = |x| + |y|.")
    print(f"Consider the points u = {u} and v = {v}.")
    print(f"The norm of u is ||u|| = |{u[0]}| + |{u[1]}| = {l1_norm(u)}")
    print(f"The norm of v is ||v|| = |{v[0]}| + |{v[1]}| = {l1_norm(v)}")
    print("Both points lie on the unit sphere.\n")

    # 2. Show that the line segment between u and v lies on the unit sphere.
    # A point on the segment is x(t) = (1-t)u + tv = (1-t, t) for t in [0,1].
    # Its L1-norm is |1-t| + |t| = (1-t) + t = 1.
    # Let's test the midpoint (t=0.5) as an example.
    t = 0.5
    midpoint = ( (1-t)*u[0] + t*v[0], (1-t)*u[1] + t*v[1] )
    print(f"The midpoint of the segment is {midpoint}.")
    print(f"The norm of the midpoint is ||{midpoint}|| = |{midpoint[0]}| + |{midpoint[1]}| = {l1_norm(midpoint)}")
    print("The entire segment lies on the unit sphere.\n")

    # 3. Calculate the length of the segment, which is ||u - v||.
    # The triangle inequality states ||u-v|| <= ||u|| + ||v|| = 1 + 1 = 2.
    # Our example achieves this maximum bound.
    u_minus_v = (u[0] - v[0], u[1] - v[1])
    length = l1_norm(u_minus_v)

    print("The length of the segment is ||u - v||.")
    print("The final equation with numbers from our example is:")
    print(f"||{u} - {v}|| = ||{u_minus_v}|| = |{u_minus_v[0]}| + |{u_minus_v[1]}| = {abs(u_minus_v[0])} + {abs(u_minus_v[1])} = {length}")

if __name__ == "__main__":
    main()