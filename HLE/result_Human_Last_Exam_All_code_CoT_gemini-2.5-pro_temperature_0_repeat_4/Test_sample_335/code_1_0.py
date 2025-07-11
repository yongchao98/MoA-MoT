import snappy
import mpmath

def solve_knot_volume():
    """
    Computes the floor of 10^6 times the simplicial volume of the knot
    K = C_{4,3}(Conway) # Wh_-^2(Eight).
    """
    # Set the precision for mpmath calculations
    mpmath.mp.dps = 50

    # --- Part 1: C_{4,3}(Conway) ---

    # The companion is the Conway knot (11n34). Get its manifold from SnapPy.
    conway_knot = snappy.Manifold('11n34')
    # Get its hyperbolic volume with high precision.
    vol_conway = mpmath.mpf(str(conway_knot.volume(bits_prec=200)))

    # v3 is the volume of the regular ideal hyperbolic tetrahedron.
    v3 = mpmath.mpf('1.0149416064096536250212025541633984385722323455512')

    # The simplicial volume of the Conway knot complement is Vol(Conway) / v3.
    norm_conway = vol_conway / v3

    # For the (4,3)-cable of the Conway knot, the winding number is 4.
    # The simplicial volume V1 is |4| * ||S^3 \ Conway||.
    w1 = 4
    V1 = abs(w1) * norm_conway

    # --- Part 2: Wh_-^2(Eight) ---

    # The companion is the figure-8 knot (4_1).
    # Its simplicial volume is exactly 2.
    norm_eight = mpmath.mpf(2)

    # For the 2-twisted negative Whitehead double, the winding number is -2.
    # The simplicial volume V2 is |-2| * ||S^3 \ Eight||.
    w2 = -2
    V2 = abs(w2) * norm_eight

    # --- Part 3: Total Volume and Final Calculation ---

    # The total simplicial volume V is the sum of V1 and V2 due to the connected sum.
    V = V1 + V2

    # The final result is floor(10^6 * V).
    result = mpmath.floor(10**6 * V)

    # --- Output the results step-by-step ---
    print("Calculation Steps:")
    print("Let V be the simplicial volume of S^3 \\ K, where K = C_{4,3}(Conway) # Wh_-^2(Eight).")
    print("V = V1 + V2, where V1 = ||S^3 \\ C_{4,3}(Conway)|| and V2 = ||S^3 \\ Wh_-^2(Eight)||.")
    
    print("\n1. Computing V1:")
    print(f"   V1 = |winding number| * ||S^3 \\ Conway||")
    print(f"   V1 = {abs(w1)} * (Vol(Conway) / v3)")
    print(f"   V1 = {abs(w1)} * ({vol_conway} / {v3})")
    print(f"   V1 = {V1}")

    print("\n2. Computing V2:")
    print(f"   V2 = |winding number| * ||S^3 \\ Eight||")
    print(f"   V2 = {abs(w2)} * {norm_eight}")
    print(f"   V2 = {V2}")

    print("\n3. Computing total volume V:")
    print(f"   V = V1 + V2 = {V1} + {V2}")
    print(f"   V = {V}")

    print("\n4. Final Calculation:")
    print(f"   The required value is floor(10^6 * V).")
    print(f"   floor(1000000 * {V}) = floor({10**6 * V})")
    print(f"   Result = {result}")

solve_knot_volume()