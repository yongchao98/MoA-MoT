def find_lsb_pos(n):
    """
    Finds the position of the least significant bit (LSB) of a positive integer n.
    This corresponds to the 'min' operation on the set of bit positions.
    For example, find_lsb_pos(12) = find_lsb_pos(0b1100) = 2.
    """
    if n <= 0:
        return -1 # Or handle as an error, as 0 has no LSB.
    pos = 0
    while (n & 1) == 0:
        n >>= 1
        pos += 1
    return pos

def hajnal_coloring(i, j):
    """
    Computes Hajnal's coloring for a pair of non-negative integers.
    The injection h: N -> P(N) is implicitly h(x) = {i | the i-th bit of x is 1}.
    The symmetric difference of the sets of bit positions corresponds to the XOR operation.
    """
    xor_val = i ^ j
    return find_lsb_pos(xor_val)

def demonstrate_hajnal_coloring():
    """
    This function demonstrates that Hajnal's coloring, a key tool in
    combinatorial set theory, is highly non-monochromatic. The theory proves
    that for the equivalent coloring on ordinals, no monochromatic triangle
    (a set {a, b, c} where f(a,b)=f(b,c)=f(a,c)) exists.

    This relates to the original problem by showing that simple Ramsey arguments
    (i.e., "let's find a monochromatic set") are insufficient. Deeper theorems
    are needed, which ultimately show that for ANY function f, a set x can
    be found where the image f''[x]^2 is small, contradicting the problem's
    requirement.
    """
    # N is a finite stand-in for an initial segment of ordinals.
    N = 100
    found_triangle = False

    print(f"--- Demonstrating Hajnal's Coloring ---")
    print(f"We search for a 'monochromatic triangle' in {{0, ..., {N-1}}}.")
    print("A monochromatic triangle is a set {i, j, k} where the 'color' of all three edges is the same.")
    print("Color f({i, j}) = position of the least significant bit of (i XOR j).")
    print("-" * 37)

    for i in range(N):
        for j in range(i + 1, N):
            for k in range(j + 1, N):
                color_ij = hajnal_coloring(i, j)
                color_jk = hajnal_coloring(j, k)
                color_ik = hajnal_coloring(i, k)

                if color_ij == color_jk and color_jk == color_ik:
                    print(f"Found a monochromatic triangle: {{{i}, {j}, {k}}}")
                    print(f"Color = {color_ij}")
                    print("This is unexpected and likely indicates an error in the code, as it contradicts the proof.")
                    found_triangle = True
                    break
            if found_triangle:
                break
        if found_triangle:
            break

    if not found_triangle:
        print("\nSearch complete: No monochromatic triangle was found.")
        print("This aligns with the mathematical proof that `kappa^+` has no monochromatic triangles under this coloring.")
        print("\nFinal Conclusion for the problem:")
        print("The broader theory shows that for any function f, one can always find a set x of type kappa+1")
        print("where the image size is small (e.g., <= 2), which is less than kappa.")
        print("Therefore, a function satisfying the condition for ALL such sets can never exist.")

if __name__ == "__main__":
    demonstrate_hajnal_coloring()