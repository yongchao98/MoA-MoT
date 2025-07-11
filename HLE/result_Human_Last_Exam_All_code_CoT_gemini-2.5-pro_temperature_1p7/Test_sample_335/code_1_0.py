import math
import snappy

def compute_simplicial_volume_value():
    """
    This function computes the value of floor(10^6 * V), where V is the
    simplicial volume of the complement of the knot K.
    K = C_{4,3}(Conway) # Wh_-^2(Eight)
    """

    # Step 1: Get the volume of a regular ideal tetrahedron, v3.
    # v3 is a fundamental constant in hyperbolic geometry.
    v3 = snappy.IdealTetrahedron().volume()

    # Step 2: Get the hyperbolic volume of the Conway knot complement.
    # The Conway knot is designated as '11n34'.
    conway_knot = snappy.Manifold("11n34")
    vol_conway = conway_knot.volume()

    # Step 3: Get the hyperbolic volume of the (-2)-twisted Whitehead link complement.
    # This is obtained by -1/2 Dehn filling on one cusp of the Whitehead link complement ('L5a1').
    # Due to the symmetry of the Whitehead link, filling either cusp gives the same result.
    whitehead_link_complement = snappy.Manifold("L5a1")
    wl_minus_2_complement = whitehead_link_complement.copy()
    # The surgery slope (-1, 2) corresponds to -1/2 filling.
    wl_minus_2_complement.dehn_fill((-1, 2), 0)
    vol_wl_minus_2 = wl_minus_2_complement.volume()

    # Step 4: Calculate the total simplicial volume V based on the formula:
    # V = (1/4 * Vol(Conway)/v3) + (Vol(WL_-2)/v3)
    # The formula uses high-precision numbers provided by snappy.
    simplicial_vol_conway = vol_conway / v3
    simplicial_vol_cabled_conway = simplicial_vol_conway / 4
    simplicial_vol_wl_minus_2 = vol_wl_minus_2 / v3
    
    V = simplicial_vol_cabled_conway + simplicial_vol_wl_minus_2
    
    # Step 5: Print the components of the calculation as requested.
    print("The final calculation is based on the formula:")
    print("V = (1/4) * (Vol_Conway / v3) + (Vol_WL_minus_2 / v3)\n")
    print(f"The values used are:")
    print(f"Vol_Conway = {vol_conway}")
    print(f"Vol_WL_minus_2 = {vol_wl_minus_2}")
    print(f"v3 = {v3}\n")
    print(f"The total simplicial volume is V = {V}\n")

    # Step 6: Compute and print the final integer result.
    result = math.floor(1000000 * V)
    print(f"The required value is floor(10^6 * V) = {result}")
    
    return result

if __name__ == "__main__":
    compute_simplicial_volume_value()