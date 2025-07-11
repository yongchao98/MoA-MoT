import collections

def solve():
    """
    Calculates the number of conjugacy classes of the extraspecial group of order 27
    and exponent 3, denoted S = 3^{1+2}_+. This number corresponds to the number
    of blocks of the group algebra kG as explained in the reasoning.
    """
    
    # S = <a, b, c | a^3=b^3=c^3=1, [a,b]=c, [a,c]=[b,c]=1>
    # Elements are represented as (i, j, k) for a^i b^j c^k, where i,j,k are in {0,1,2}.
    
    # Multiplication rule: (i1,j1,k1) * (i2,j2,k2) = (i1+i2, j1+j2, k1+k2+j1*i2) mod 3
    def multiply(g1, g2):
        i1, j1, k1 = g1
        i2, j2, k2 = g2
        return ((i1 + i2) % 3, (j1 + j2) % 3, (k1 + k2 + j1 * i2) % 3)

    # Inverse rule: inv(i,j,k) = (-i, -j, -k+i*j) mod 3
    def inverse(g):
        i, j, k = g
        return ((-i) % 3, (-j) % 3, (-k + i * j) % 3)

    # Generate all 27 elements of the group S
    elements = []
    for i in range(3):
        for j in range(3):
            for k in range(3):
                elements.append((i, j, k))

    # Algorithm to find conjugacy classes
    unclassified = set(elements)
    classes = []
    while unclassified:
        x = unclassified.pop()
        current_class = set()
        for g in elements:
            # conjugate(g, x) = g * x * inv(g)
            conj = multiply(multiply(g, x), inverse(g))
            current_class.add(conj)
        classes.append(frozenset(current_class))
        unclassified -= current_class

    # Print the results
    print("The number of blocks of kG is equal to the number of conjugacy classes of S.")
    print("Calculating the conjugacy classes of S = 3^{1+2}_+:")
    print("-" * 30)
    
    # Sort classes for consistent output, primarily by size then by elements
    sorted_classes = sorted([sorted(list(c)) for c in classes])

    class_counts = collections.defaultdict(int)
    for i, cls in enumerate(sorted_classes):
        size = len(cls)
        class_counts[size] += 1
        print(f"Class {i+1} (size {size}):")
        # To avoid overly long output, print at most 5 elements
        elements_to_print = ", ".join(map(str, cls[:5]))
        if len(cls) > 5:
            elements_to_print += ", ..."
        print(f"  {{{elements_to_print}}}")
    
    print("-" * 30)
    print("Summary of class sizes:")
    
    total_classes = 0
    equation_parts = []
    for size, count in sorted(class_counts.items()):
        print(f"- {count} classes of size {size}")
        total_classes += count
        equation_parts.append(str(count))

    # As per the prompt, outputting the numbers in the final equation.
    # Here the "equation" is the sum of the counts of classes of each size.
    # The calculation is trivial in this case as there are only two sizes.
    count_size_1 = class_counts.get(1, 0)
    count_size_3 = class_counts.get(3, 0)
    print(f"\nThe total number of conjugacy classes is the sum of counts for each size.")
    print(f"Number of classes = (Number of classes of size 1) + (Number of classes of size 3)")
    print(f"Final Equation: {count_size_1} + {count_size_3} = {total_classes}")
    print(f"\nTherefore, the number of blocks of kG is {total_classes}.")

solve()