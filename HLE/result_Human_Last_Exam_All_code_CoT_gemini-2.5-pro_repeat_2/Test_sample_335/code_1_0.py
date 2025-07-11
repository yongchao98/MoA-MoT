import math

def compute_knot_volume():
    """
    Computes the floor of 10^6 * V, where V is the simplicial volume of the
    complement of the knot K := C_{4,3}(Conway) # Wh_-^2(Eight).
    """
    
    # The knot K is a connected sum K1 # K2.
    # K1 = C_{4,3}(Conway)
    # K2 = Wh_-^2(Eight)
    # The simplicial volume V is additive over connected sums: V = V(K1) + V(K2).
    
    print("Step 1: Calculate the simplicial volume of the complement of K1 = C_{4,3}(Conway).")
    # For a cable knot C_{p,q}(J), the volume is |p| * ||S^3 \ J||.
    p_cable = 4
    # The Conway knot is a slice knot, so the simplicial volume of its complement is 0.
    V_conway = 0
    V_K1 = p_cable * V_conway
    print(f"The formula for a cable knot is ||S^3 \\ C_{{p,q}}(J)|| = |p| * ||S^3 \\ J||.")
    print(f"For K1, p = {p_cable} and J is the Conway knot.")
    print(f"The simplicial volume of the Conway knot complement is {V_conway} because it is a slice knot.")
    print(f"So, ||S^3 \\ K1|| = {p_cable} * {V_conway} = {V_K1}")
    print("-" * 20)

    print("Step 2: Calculate the simplicial volume of the complement of K2 = Wh_-^2(Eight).")
    # For a Whitehead double Wh_n(J), the volume is ||S^3 \ J||.
    # The companion knot J is the figure-8 knot.
    # The simplicial volume of the figure-8 knot complement is 2.
    V_eight = 2
    V_K2 = V_eight
    print("The formula for a Whitehead double is ||S^3 \\ Wh_n(J)|| = ||S^3 \\ J||.")
    print("For K2, J is the figure-8 knot.")
    print(f"The simplicial volume of the figure-8 knot complement is a known value: {V_eight}.")
    print(f"So, ||S^3 \\ K2|| = {V_K2}")
    print("-" * 20)

    print("Step 3: Calculate the total simplicial volume V.")
    # Total simplicial volume V = V(K1) + V(K2)
    V = V_K1 + V_K2
    print(f"V = ||S^3 \\ K1|| + ||S^3 \\ K2||")
    print(f"V = {V_K1} + {V_K2} = {V}")
    print("-" * 20)

    print("Step 4: Compute the final value floor(10^6 * V).")
    factor = 10**6
    final_result = math.floor(factor * V)
    print(f"The required value is floor({factor} * V).")
    print(f"floor({factor} * {V}) = floor({factor * V}) = {final_result}")
    
if __name__ == '__main__':
    compute_knot_volume()
