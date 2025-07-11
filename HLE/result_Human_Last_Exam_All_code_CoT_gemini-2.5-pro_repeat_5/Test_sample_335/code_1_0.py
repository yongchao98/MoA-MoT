import snappy
import mpmath

def solve_knot_volume():
    """
    This function calculates the floor of 10^6 * V, where V is the simplicial volume
    of the complement of the knot K = C_{4,3}(Conway) # Wh_-^2(Eight).
    """
    # Set high precision for calculations
    mpmath.mp.dps = 50

    # The volume of a regular ideal tetrahedron
    # v3 = Im(Li_2(e^{i*pi/3}))
    v3 = mpmath.imag(mpmath.li(2, mpmath.exp(mpmath.j * mpmath.pi / 3)))

    # Part 1: C_{4,3}(Conway)
    # The simplicial volume of the Conway knot complement is the simplicial volume
    # of the P(-2,3,7) pretzel knot complement.
    # SnapPy can identify this knot and compute its volume.
    knot_p_neg2_3_7 = snappy.Manifold('p(-2,3,7)')
    vol_conway_pattern = knot_p_neg2_3_7.volume()
    simplicial_vol_conway = vol_conway_pattern / v3

    # The simplicial volume of the (4,3)-cable is 4 times that of the companion.
    simplicial_vol_k1 = 4 * simplicial_vol_conway

    # Part 2: Wh_-^2(Eight)
    # The simplicial volume of a Whitehead double is independent of the companion knot.
    # We compute the volume for the double of the unknot.
    # This corresponds to (-1/2)-Dehn surgery on a component of the Whitehead link (L5a1).
    whitehead_link_exterior = snappy.Manifold('L5a1')
    # Perform (-1/2)-surgery on the second component (index 1)
    knot_k2_exterior = whitehead_link_exterior.dehn_fill((-1, 2), 1)
    vol_k2 = knot_k2_exterior.volume()
    simplicial_vol_k2 = vol_k2 / v3

    # Total simplicial volume V
    V = simplicial_vol_k1 + simplicial_vol_k2
    
    # Print the detailed calculation
    print("Step 1: Calculate the simplicial volume of the C_{4,3}(Conway) part.")
    print(f"Volume of P(-2,3,7) complement: Vol_p = {vol_conway_pattern}")
    print(f"Volume of ideal tetrahedron: v3 = {v3}")
    print(f"Simplicial volume of Conway complement: ||S^3\\Conway|| = Vol_p / v3 = {simplicial_vol_conway}")
    print(f"Simplicial volume of C_{4,3}(Conway): V1 = 4 * ||S^3\\Conway|| = {simplicial_vol_k1}\n")

    print("Step 2: Calculate the simplicial volume of the Wh_-^2(Eight) part.")
    print(f"Volume of Wh_-^2(unknot) complement: Vol_k2 = {vol_k2}")
    print(f"Simplicial volume of Wh_-^2(Eight): V2 = Vol_k2 / v3 = {simplicial_vol_k2}\n")

    print("Step 3: Calculate the total simplicial volume V.")
    print(f"V = V1 + V2 = {simplicial_vol_k1} + {simplicial_vol_k2} = {V}\n")
    
    # Final computation
    final_value = mpmath.floor(10**6 * V)
    
    print("Step 4: Calculate the final answer.")
    print(f"The total simplicial volume is V = ({4} * {vol_conway_pattern} + {vol_k2}) / {v3}")
    print(f"V approx {V}")
    print(f"The required value is floor(10^6 * V) = floor({10**6 * V})")
    print(f"Final answer: {final_value}")

if __name__ == '__main__':
    solve_knot_volume()
