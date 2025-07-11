def solve():
    """
    This function calculates and formats the properties for the three given CFGs.
    """

    # Properties for X_1: Hilb^11(A^3)
    # Type: Scheme
    # Separated: Yes
    # Universally closed: No
    # Irreducible: No
    # Dimension: 11 * 3 = 33
    x1_props = ['S', 's', 11 * 3]
    x1_str = f"[{','.join(map(str, x1_props))}]"

    # Properties for X_2: [(\A^4 \ V(xy-zw))/C*]
    # Type: DM stack
    # Separated: Yes
    # Universally closed: No
    # Irreducible: Yes
    # Dimension: 4 - 1 = 3
    x2_props = ['DM', 's', 'irr', 4 - 1]
    x2_str = f"[{','.join(map(str, x2_props))}]"

    # Properties for X_3: Pic(C_0) for a genus 7 curve
    # Type: Algebraic stack
    # Separated: Yes
    # Universally closed: No
    # Irreducible: No
    # Dimension: genus = 7
    g = 7
    x3_props = ['A', 's', g]
    x3_str = f"[{','.join(map(str, x3_props))}]"

    # Combine and print the results
    final_answer = ' '.join([x1_str, x2_str, x3_str])
    print(final_answer)

solve()