import sys

def calculate_riemann_components_kahler():
    """
    Calculates the number of independent components of the Riemann tensor
    on a Kähler manifold of a given complex dimension.
    """
    try:
        m_str = input("Enter the complex dimension of the Kähler manifold (m): ")
        m = int(m_str)
        if m < 1:
            print("Error: The complex dimension must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Please enter a valid integer for the dimension.", file=sys.stderr)
        return

    # Step 1: The space of Kähler curvature tensors at a point is isomorphic
    # to the space of symmetric bilinear forms on the Lie Algebra u(m), S^2(u(m)).
    # We first need the dimension of u(m).
    # dim(u(m)) = m^2
    d = m**2

    # Step 2: The dimension of the space of symmetric bilinear forms on a
    # vector space of dimension d is d*(d+1)/2.
    num_components = (d * (d + 1)) // 2

    print("\n--- Calculation Steps ---")
    print(f"For a Kähler manifold of complex dimension m = {m}:")
    print(f"1. The relevant Lie algebra is u(m), which has a real dimension d = m^2.")
    print(f"   d = {m}^2 = {d}")
    print("\n2. The number of independent components of the Riemann tensor is the dimension of the space of")
    print(f"   symmetric bilinear forms on u(m), which is given by the formula: d*(d+1)/2.")
    print(f"   Substituting d = {d}, we get: ({d} * ({d} + 1)) / 2 = {num_components}")

    print("\n--- Final Equation ---")
    d_plus_1 = d + 1
    # Print each number in the final equation as requested.
    print(f"({m**2} * ({m**2} + 1)) / 2 = ({d} * {d_plus_1}) / 2 = {num_components}")
    print(f"\nThe number of independent real entries for the Riemann tensor is: {num_components}")


if __name__ == '__main__':
    calculate_riemann_components_kahler()
