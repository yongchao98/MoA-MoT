def solve():
    """
    This function provides the properties of the three given categories fibered in groupoids.
    """
    
    # Properties for X1: The Hilbert scheme Hilb_11(A^3)
    # Type: Scheme (S), separated (s), dim = 11 * 3 = 33
    x1_props = "[S,s,33]"

    # Properties for X2: The quotient stack [(A^4 \ V(xy-zw))/C*]
    # Type: DM stack (DM), separated (s), irreducible (irr), dim = 4 - 1 = 3
    x2_props = "[DM,s,irr,3]"
    
    # Properties for X3: The Picard stack Pic(C_0) for a genus 7 curve
    # Type: Algebraic stack (A), separated (s), dim = g - 1 = 7 - 1 = 6
    x3_props = "[A,s,6]"
    
    final_answer = " ".join([x1_props, x2_props, x3_props])
    print(final_answer)

solve()