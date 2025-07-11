import numpy as np

def identify_manifold_from_heegaard():
    """
    Identifies the 3-manifold by calculating the order of its first homology group
    from a known presentation of its fundamental group.

    The Heegaard diagram represents the Poincaré homology sphere. Its fundamental group
    can be presented as <s, t | (st)^2 = s^3, s^3 = t^5>.

    The first homology group H1(M) is the abelianization of the fundamental group.
    The relations in additive notation are:
    -s + 2t = 0
    3s - 5t = 0
    """

    # The relation matrix is formed by the coefficients of the abelianized relations.
    M = np.array([
        [-1, 2],
        [3, -5]
    ])

    print("The manifold represented by the Heegaard diagram is the Poincaré homology sphere.")
    print("A known presentation for its fundamental group is <s, t | (st)^2 = s^3, s^3 = t^5>.")
    print("The abelianized relations yield the following integer matrix:")
    print(M)
    print("\nTo find the order of the first homology group H1(M), we calculate the absolute value of the determinant of this matrix.")
    
    # Calculate determinant using integer arithmetic to be precise
    a, b = M[0, 0], M[0, 1]
    c, d = M[1, 0], M[1, 1]
    det_M = a * d - b * c

    print(f"The determinant is calculated as: det(M) = ({a})*({d}) - ({b})*({c})")
    print(f"det(M) = {det_M}")

    order_H1 = abs(det_M)

    print(f"\nThe order of the first homology group is |det(M)| = {order_H1}.")

    if order_H1 == 1:
        print("Since the order is 1, the group is trivial. This confirms that the manifold is a homology sphere, as expected.")
    else:
        print(f"The homology group has order {order_H1}, so this manifold is not a homology sphere.")

if __name__ == '__main__':
    identify_manifold_from_heegaard()