import collections

def group_multiply(g1, g2):
    """
    Computes the product of two elements in the extraspecial group S = 3^{1+2}_+.
    Elements are represented as tuples (i, j, k) for a^i b^j c^k.
    The multiplication rule is derived from ba = abc^2.
    (a^i b^j c^k) * (a^x b^y c^z) = a^(i+x) b^(j+y) c^(k+z+2jx) mod 3.
    """
    i, j, k = g1
    x, y, z = g2
    # The calculation is done in the exponents, modulo 3.
    new_i = (i + x) % 3
    new_j = (j + y) % 3
    new_k = (k + z + 2 * j * x) % 3
    return (new_i, new_j, new_k)

def group_inverse(g):
    """
    Computes the inverse of an element in S.
    The inverse of a^i b^j c^k is a^(-i) b^(-j) c^(-k+2ij) mod 3.
    """
    i, j, k = g
    # Exponents are modulo 3.
    inv_i = (-i) % 3
    inv_j = (-j) % 3
    inv_k = (-k + 2 * i * j) % 3
    return (inv_i, inv_j, inv_k)

def main():
    """
    Main function to compute the number of conjugacy classes for S.
    """
    # Generate all 27 elements of the group S.
    all_elements = set()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                all_elements.add((i, j, k))

    remaining_elements = all_elements.copy()
    class_sizes = []

    # Partition the group into conjugacy classes.
    while remaining_elements:
        # Pick an element to start a new class.
        g = remaining_elements.pop()

        # Compute its conjugacy class C(g) = {hgh^-1 for h in S}.
        conjugacy_class = set()
        for h in all_elements:
            h_inv = group_inverse(h)
            h_g = group_multiply(h, g)
            conjugate = group_multiply(h_g, h_inv)
            conjugacy_class.add(conjugate)

        # Store the size of the class.
        class_sizes.append(len(conjugacy_class))
        
        # Remove elements of this class from the remaining set.
        # Note: 'g' was already popped.
        remaining_elements.difference_update(conjugacy_class)
    
    # Sort class sizes for a clean presentation.
    class_sizes.sort()

    # The number of blocks is the sum of the class sizes, which is the number of classes.
    num_blocks = len(class_sizes)
    
    # Format the sizes into an equation.
    equation_str = " + ".join(map(str, class_sizes))

    print(f"The number of blocks is determined by the number of conjugacy classes of S.")
    print(f"The sizes of the conjugacy classes are: {class_sizes}")
    print(f"The total number of classes is the sum of counts for each class size.")
    print(f"Thus, the number of blocks of kG is:")
    print(f"{equation_str} = {num_blocks}")


if __name__ == "__main__":
    main()