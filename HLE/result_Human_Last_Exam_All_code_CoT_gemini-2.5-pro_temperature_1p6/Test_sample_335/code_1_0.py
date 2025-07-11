import math

def compute_knot_volume():
    """
    Computes the floor of 10^6 times the simplicial volume of the specified knot.

    The knot is K = C_{4,3}(Conway) # Wh_-^2(Eight).
    Its simplicial volume V = ||S^3 \ K|| is calculated as:
    V = ||S^3 \ C_{4,3}(Conway)|| + ||S^3 \ Wh_-^2(Eight)||  (additivity over connected sum #)
      = (4 * ||S^3 \ Conway||) + 0                            (cabling and Whitehead double properties)
      = 4 * Vol_JSJ(S^3 \ Conway) / v3                         (relation between simplicial and hyperbolic volume)
    """

    # The JSJ-volume of the Conway knot complement is the hyperbolic volume of its mutant, 11n42.
    vol_conway_jsj = 2.828122088316131

    # v3 is the volume of the regular ideal tetrahedron.
    v3 = 1.0149416064096536

    # p is the winding number from the (4,3)-cable pattern.
    p = 4

    # Calculate the simplicial volume V
    V = p * vol_conway_jsj / v3

    print("The formula for the simplicial volume V is: p * (Vol_JSJ(Conway) / v3)")
    print(f"where p (winding number) = {p}")
    print(f"Vol_JSJ(Conway) (hyperbolic volume of its mutant 11n42) = {vol_conway_jsj}")
    print(f"v3 (volume of regular ideal tetrahedron) = {v3}")
    print(f"V = {p} * ({vol_conway_jsj} / {v3})")
    print(f"Calculated simplicial volume V = {V}")

    # Compute the final result
    result = math.floor(10**6 * V)

    print(f"\nThe value to compute is floor(10^6 * V)")
    print(f"10^6 * V = {10**6 * V}")
    print(f"The final result is: {result}")

compute_knot_volume()