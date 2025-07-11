import math

def solve():
    """
    Calculates the maximum number of grid squares a triangle with sides
    18, 18, and 18*sqrt(2) can pass through without its perimeter
    containing any lattice points.

    The optimal orientation is placing the hypotenuse parallel to an axis.
    """

    # Side lengths
    leg_length = 18
    hypotenuse_length = 18 * math.sqrt(2)
    half_hypotenuse_projection = 9 * math.sqrt(2)

    # Number of grid lines crossed by each segment in the optimal orientation.
    # We use floor() because a line from x=0.1 to x=25.5 crosses lines x=1, 2, ..., 25.

    # 1. Squares crossed by the hypotenuse (Segment AB)
    # Placed nearly parallel to the x-axis.
    # Crosses vertical lines only.
    v_cross_AB = math.floor(hypotenuse_length)
    h_cross_AB = 0
    squares_AB = 1 + v_cross_AB + h_cross_AB

    # 2. Squares crossed by the first leg (Segment AC)
    # This leg connects a hypotenuse endpoint to the right-angle vertex.
    # Its projections on both axes are 9*sqrt(2).
    v_cross_AC = math.floor(half_hypotenuse_projection)
    h_cross_AC = math.floor(half_hypotenuse_projection)
    squares_AC = 1 + v_cross_AC + h_cross_AC

    # 3. Squares crossed by the second leg (Segment BC)
    # Projections are also 9*sqrt(2), but over a different grid interval.
    # Vertical crossings: from x=ceil(9*sqrt(2)) to floor(18*sqrt(2))
    # Horizontal crossings: from y=0 to floor(9*sqrt(2))
    v_cross_BC = math.floor(hypotenuse_length) - math.floor(half_hypotenuse_projection)
    h_cross_BC = math.floor(half_hypotenuse_projection)
    squares_BC = 1 + v_cross_BC + h_cross_BC

    # 4. Overlaps between segments
    # The squares containing the vertices are common to the segments meeting there.
    # Overlap AB-AC: The square containing vertex A.
    overlap_AB_AC = 1
    # Overlap AC-BC: The square containing vertex C.
    overlap_AC_BC = 1
    # Overlap AB-BC: Near vertex B. The leg BC enters the y=0 strip that AB lies in,
    # crossing two squares in that strip before ending.
    overlap_AB_BC = 2

    # 5. Total unique squares using Principle of Inclusion-Exclusion
    total_k = squares_AB + squares_AC + squares_BC - (overlap_AB_AC + overlap_AC_BC + overlap_AB_BC)

    print("Calculation Steps:")
    print(f"1. Squares for hypotenuse (AB): {squares_AB}")
    print(f"2. Squares for first leg (AC): {squares_AC}")
    print(f"3. Squares for second leg (BC): {squares_BC}")
    print(f"4. Overlapping squares: {overlap_AB_AC} (at A), {overlap_AC_BC} (at C), {overlap_AB_BC} (at B)")
    print("\nFinal Equation:")
    print(f"k = {squares_AB} + {squares_AC} + {squares_BC} - ({overlap_AB_AC} + {overlap_AC_BC} + {overlap_AB_BC})")
    print(f"k = {total_k}")

solve()
<<<73>>>