import snappy
import math

def compute_simplicial_volume_floor():
    """
    Computes the value of floor(10^6 * V), where V is the simplicial volume
    of the complement of the knot K = C(4,3)(Conway) # Wh_-^2(Eight).
    """
    # Set high precision for numerical calculations to ensure accuracy.
    snappy.Manifold.set_precision(bits=100)

    print("The final equation for the simplicial volume V is:")
    print("V = (Vol_Conway + Vol_Eight + Vol_Pattern) / v_3")
    print("\n--- Calculating the components of the equation ---")

    # 1. Calculate v_3, the volume of a regular ideal hyperbolic tetrahedron.
    # This is half the volume of the figure-eight knot complement.
    eight_knot_manifold = snappy.Manifold('4_1')
    vol_eight = eight_knot_manifold.volume()
    v_3 = vol_eight / 2
    print(f"Vol_Eight (Volume of S^3 \\ Figure-8 knot): {vol_eight}")
    print(f"v_3 (Volume of regular ideal tetrahedron): {v_3}")

    # 2. Calculate the volume of the Conway knot complement.
    # The Conway knot is denoted as 11n34 in the Rolfsen table.
    conway_knot_manifold = snappy.Manifold('11n34')
    vol_conway = conway_knot_manifold.volume()
    print(f"Vol_Conway (Volume of S^3 \\ Conway knot): {vol_conway}")

    # 3. Calculate the volume of the Whitehead pattern space.
    # This space is obtained by -1/2 Dehn surgery on one component of the
    # negative Whitehead link complement (L5a1).
    whitehead_link_manifold = snappy.Manifold('L5a1')
    # The two components of the Whitehead link are symmetric, so we can fill either cusp.
    pattern_space_manifold = whitehead_link_manifold.copy()
    pattern_space_manifold.dehn_fill(slopes=(-1, 2), cusp_index=0)
    vol_pattern = pattern_space_manifold.volume()
    print(f"Vol_Pattern (Volume of the Wh_-^2 pattern space): {vol_pattern}")

    # 4. Sum the hyperbolic volumes of all hyperbolic JSJ pieces.
    total_hyperbolic_volume = vol_conway + vol_eight + vol_pattern
    print("\n--- Final Calculation ---")
    print(f"Total hyperbolic volume = {vol_conway} + {vol_eight} + {vol_pattern}")
    print(f"                        = {total_hyperbolic_volume}")

    # 5. Compute the total simplicial volume V.
    simplicial_volume_V = total_hyperbolic_volume / v_3
    print(f"Simplicial Volume V = {total_hyperbolic_volume} / {v_3}")
    print(f"                    = {simplicial_volume_V}")

    # 6. Compute the final result as requested.
    final_answer = math.floor(10**6 * simplicial_volume_V)
    print(f"\nValue to compute: floor(10^6 * V) = floor(1000000 * {simplicial_volume_V})")
    print(f"Final result: {final_answer}")
    
    return final_answer

if __name__ == '__main__':
    result = compute_simplicial_volume_floor()
    print(f"\n<<<{result}>>>")
