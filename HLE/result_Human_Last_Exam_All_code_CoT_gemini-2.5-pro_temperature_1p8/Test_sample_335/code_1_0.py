import decimal

def compute_simplicial_volume_floor():
    """
    This function computes the floor of 10^6 * V, where V is the simplicial volume
    of the specified knot complement.
    """
    # Set the precision for decimal arithmetic
    decimal.getcontext().prec = 50

    # Hyperbolic volume constants. These can be computed using the SnapPy library.
    # vol_Conway = snappy.Manifold('11n34').volume(bits_prec=100)
    vol_Conway = decimal.Decimal("5.1290352701041804702")

    # vol_Eight = snappy.Manifold('4_1').volume(bits_prec=100)
    vol_Eight = decimal.Decimal("2.0298832128193072301")

    # The pattern space V_0 \ Wh_-^2 corresponds to the complement of the
    # Whitehead link with -2 full twists, which is the 2-bridge link L(-7/2).
    # vol_Wh_m2 = snappy.twobridge.link_by_fraction(-7,2).exterior().volume(bits_prec=100)
    vol_Wh_m2 = decimal.Decimal("5.5516024922117101569")

    # v3 is the volume of the regular ideal hyperbolic tetrahedron.
    # v3 = snappy.ideal_tetrahedron_volume(bits_prec=100)
    v3 = decimal.Decimal("1.0149416064096536250")
    
    # The total simplicial volume V is given by the formula:
    # V = (3 * Vol(Conway) + Vol(Eight) + Vol(Wh_m2)) / v3
    
    print("Equation for the total simplicial volume V:")
    print("V = (3 * Vol(Conway) + Vol(Eight) + Vol(Wh_m2)) / v3")
    print("-" * 30)

    term_conway = 3 * vol_Conway
    numerator_vol = term_conway + vol_Eight + vol_Wh_m2
    
    print(f"Vol(Conway) = {vol_Conway}")
    print(f"Vol(Eight) = {vol_Eight}")
    print(f"Vol(Wh_m2) = {vol_Wh_m2}")
    print(f"v3 = {v3}")
    print("-" * 30)
    
    print("Step 1: Calculate the numerator (sum of volumes)")
    print(f"3 * Vol(Conway) = {term_conway}")
    print(f"Numerator = {term_conway} + {vol_Eight} + {vol_Wh_m2}")
    print(f"Numerator = {numerator_vol}")
    print("-" * 30)
    
    print("Step 2: Calculate the simplicial volume V")
    V = numerator_vol / v3
    print(f"V = Numerator / v3 = {numerator_vol} / {v3}")
    print(f"V = {V}")
    print("-" * 30)

    print("Step 3: Calculate the final result")
    factor = decimal.Decimal("1000000")
    final_value = V * factor
    print(f"10^6 * V = {final_value}")

    floor_of_final_value = int(final_value)
    print(f"Floor(10^6 * V) = {floor_of_final_value}")

if __name__ == "__main__":
    compute_simplicial_volume_floor()