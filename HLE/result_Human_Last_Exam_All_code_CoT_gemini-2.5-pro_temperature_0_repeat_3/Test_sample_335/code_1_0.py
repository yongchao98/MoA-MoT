import snappy
import math

def solve_knot_volume():
    """
    Computes the simplicial volume of the specified knot complement and the final requested value.
    """
    # Set precision for snappy computations to ensure accuracy.
    snappy.Manifold.set_precision(200)

    # Step 1: Obtain the required hyperbolic volumes using snappy.

    # The figure-8 knot is '4_1'. Its volume is 2*v_3.
    eight_knot_complement = snappy.Manifold('4_1')
    vol_eight = eight_knot_complement.volume()
    
    # v_3 is the volume of the regular ideal hyperbolic tetrahedron.
    v_3 = vol_eight / 2

    # The Conway knot is '11n34' in the Rolfsen table.
    conway_knot_complement = snappy.Manifold('11n34')
    vol_conway = conway_knot_complement.volume()

    # The pattern space for the 2-twisted negative Whitehead double is obtained by
    # Dehn filling on one component of the Whitehead link ('L5a1').
    # The notation Wh_-^2 corresponds to n=-2 twists, which implies a (-2,1)-Dehn surgery.
    whitehead_link_complement = snappy.Manifold('L5a1')
    wh_pattern_space = whitehead_link_complement.copy()
    wh_pattern_space.dehn_fill((-2, 1), 0)
    vol_wh_pattern = wh_pattern_space.volume()

    # Step 2: Calculate the total simplicial volume V using the formula.
    # V = (Vol(Conway) + Vol(Eight) + Vol(M_Wh)) / v_3
    total_hyperbolic_volume = vol_conway + vol_eight + vol_wh_pattern
    V = total_hyperbolic_volume / v_3

    # Step 3: Compute the final result, floor(10^6 * V).
    final_answer = math.floor(10**6 * V)

    # Step 4: Print the components of the final equation and the result.
    print("The final equation for the simplicial volume V is:")
    print(f"V = (Vol(Conway) + Vol(Eight) + Vol(M_Wh)) / v_3")
    print("\nWhere the numerical values are:")
    print(f"Vol(Conway) = {vol_conway}")
    print(f"Vol(Eight) = {vol_eight}")
    print(f"Vol(M_Wh) = {vol_wh_pattern}")
    print(f"v_3 = {v_3}")
    
    print("\nSubstituting the values:")
    print(f"V = ({vol_conway} + {vol_eight} + {vol_wh_pattern}) / {v_3}")
    print(f"V = {total_hyperbolic_volume} / {v_3}")
    print(f"V = {V}")

    print("\nThe final computation is:")
    print(f"floor(10^6 * V) = floor({10**6 * V}) = {final_answer}")

solve_knot_volume()