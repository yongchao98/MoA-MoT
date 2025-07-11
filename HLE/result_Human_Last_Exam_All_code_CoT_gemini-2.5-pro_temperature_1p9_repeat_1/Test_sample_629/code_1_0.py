import numpy as np
import itertools

def perm_multiply(p1, p2):
    """Multiplies two permutations p1 and p2 (p1 after p2)."""
    n = len(p2)
    result = [0] * n
    for i in range(n):
        result[i] = p1[p2[i]]
    return tuple(result)

def get_matrix(sigma_x, sigma_o):
    """Converts a diagram (sigma_x, sigma_o) to a matrix representation."""
    n = len(sigma_x)
    M = np.zeros((n, n), dtype=int)
    for j in range(n):
        M[sigma_x[j], j] = 1  # 'X' marker
        M[sigma_o[j], j] = -1 # 'O' marker
    return M

def find_unique_diagrams():
    """
    Calculates the number of unique minimal grid diagrams for the left-hand trefoil
    knot up to rotation.
    """
    n = 3
    # The cyclic permutation for the left-hand trefoil
    lh_trefoil_cycle = (2, 0, 1)

    # 1. Generate all 6 permutations for a 3x3 grid.
    base_perms = list(itertools.permutations(range(n)))
    
    # 2. Generate the 6 grid diagrams for the left-hand trefoil.
    lh_diagram_perms = []
    for p in base_perms:
        sigma_x = p
        sigma_o = perm_multiply(lh_trefoil_cycle, sigma_x)
        lh_diagram_perms.append((sigma_x, sigma_o))

    # 3. Convert diagrams to matrices.
    matrices = [get_matrix(sx, so) for sx, so in lh_diagram_perms]
    # We assign a name to each for clarity
    named_matrices = {f"D{i+1}": m for i, m in enumerate(matrices)}
    
    print(f"The minimal grid number for the trefoil knot is {n}.")
    print(f"Found {len(matrices)} total grid diagrams for the left-hand trefoil knot.")
    print("\nChecking for rotational equivalence...")
    
    # 4. Group matrices into rotational equivalence classes.
    unclassified_indices = list(range(len(matrices)))
    equivalence_classes = []
    
    while unclassified_indices:
        base_idx = unclassified_indices.pop(0)
        base_matrix = matrices[base_idx]
        new_class_indices = {base_idx}
        
        # Generate the 4 rotations of the base matrix
        rotations = [np.rot90(base_matrix, k=k) for k in range(4)]
        
        # Find other matrices that are rotationally equivalent
        for i in list(unclassified_indices):
            for rot_matrix in rotations:
                if np.array_equal(matrices[i], rot_matrix):
                    new_class_indices.add(i)
                    unclassified_indices.remove(i)
                    break
        equivalence_classes.append(new_class_indices)

    # 5. Print the results
    print("\nThe diagrams partition into the following classes:")
    for i, class_indices in enumerate(equivalence_classes):
        class_names = [f"D{j+1}" for j in sorted(list(class_indices))]
        print(f"Equivalence Class {i+1}: {', '.join(class_names)}")
    
    num_unique = len(equivalence_classes)
    
    print(f"\nFinal count:")
    print("The final count is the number of equivalence classes.")
    # The problem asks to output the numbers in the final equation.
    # In this case, the total diagrams are partitioned into classes.
    # We can represent this as the sum of the sizes of the classes equaling the total.
    class_sizes = [len(c) for c in equivalence_classes]
    equation_parts = [str(s) for s in class_sizes]
    print(f"{' + '.join(equation_parts)} = {sum(class_sizes)}")
    print(f"The number of unique grid diagrams is {num_unique}.")


if __name__ == '__main__':
    find_unique_diagrams()