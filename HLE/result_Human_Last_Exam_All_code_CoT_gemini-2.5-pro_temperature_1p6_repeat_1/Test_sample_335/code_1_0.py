import snappy
import math

# This script calculates the simplicial volume of the complement of the knot
# K = C_{4,3}(Conway) # Wh_-^2(Eight).
# The calculation relies on properties of the simplicial volume and hyperbolic volumes
# of specific knot complements, which are obtained using the SnapPy library.

try:
    # Step 1: Obtain the necessary hyperbolic volumes.

    # v3 is the volume of the regular ideal hyperbolic tetrahedron. It is a fundamental
    # constant in the theory of hyperbolic 3-manifolds.
    v3 = snappy.IdealTetrahedron.volume()

    # The Conway knot is denoted as '11n34' in the Rolfsen table and SnapPy.
    M_conway = snappy.Manifold('11n34')
    vol_conway = M_conway.volume()

    # The figure-8 knot is denoted as '4_1'.
    M_eight = snappy.Manifold('4_1')
    vol_eight = M_eight.volume()

    # The 2-twisted negative Whitehead double of the unknot, Wh_-^2(Unknot),
    # is topologically the knot '5_2'.
    M_5_2 = snappy.Manifold('5_2')
    vol_5_2 = M_5_2.volume()

    # Step 2: Calculate the total simplicial volume V.
    # The total simplicial volume is the sum of the hyperbolic volumes of the
    # constituent hyperbolic pieces, all divided by v3.
    total_hyperbolic_volume = vol_conway + vol_eight + vol_5_2
    V = total_hyperbolic_volume / v3

    # Step 3: Print the breakdown of the calculation as requested.
    print("The simplicial volume V is calculated using the formula:")
    print("V = (Vol(S^3 \\ Conway) + Vol(S^3 \\ Eight) + Vol(S^3 \\ 5_2)) / v3\n")

    print("The numerical values for the components are:")
    print(f"Volume of Conway knot complement: Vol(S^3 \\ Conway) = {vol_conway}")
    print(f"Volume of Figure-8 knot complement: Vol(S^3 \\ Eight) = {vol_eight}")
    print(f"Volume of 5_2 knot complement: Vol(S^3 \\ 5_2) = {vol_5_2}")
    print(f"Volume of regular ideal tetrahedron: v3 = {v3}\n")

    print("Substituting these values into the formula:")
    print(f"V = ({vol_conway} + {vol_eight} + {vol_5_2}) / {v3}")
    print(f"V = {total_hyperbolic_volume} / {v3}")
    print(f"V = {V}\n")

    # Step 4: Compute the final requested value, floor(10^6 * V).
    result = math.floor(10**6 * V)

    print(f"The value of 10^6 * V is: {10**6 * V}")
    print(f"The final result, floor(10^6 * V), is: {result}")

except ImportError:
    print("This script requires the 'snappy' library.")
    print("Please install it, for example using: pip install snappy")
except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure your snappy installation is working correctly.")
