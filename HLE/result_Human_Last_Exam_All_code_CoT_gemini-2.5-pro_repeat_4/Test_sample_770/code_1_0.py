import itertools

def get_permutation_sign(p):
    """
    Calculates the sign of a permutation.
    The sign is (-1)^(n - number of cycles).
    A permutation is even if the sign is +1.
    """
    n = len(p)
    visited = [False] * n
    num_cycles = 0
    for i in range(n):
        if not visited[i]:
            num_cycles += 1
            j = i
            while not visited[j]:
                visited[j] = True
                j = p[j]
    return (-1)**(n - num_cycles)

def multiply_perms(p1, p2):
    """Computes the composition of two permutations, p1 o p2."""
    return tuple(p1[i] for i in p2)

def invert_perm(p):
    """Computes the inverse of a permutation."""
    n = len(p)
    inv = [0] * n
    for i in range(n):
        inv[p[i]] = i
    return tuple(inv)

def calculate_conjugacy_classes(group):
    """
    Partitions a group into its conjugacy classes.
    Returns a list of frozensets, where each frozenset is a conjugacy class.
    """
    if not group:
        return []
    
    group_elements = set(group)
    conjugacy_classes = []
    
    while group_elements:
        g = group_elements.pop()
        
        # Compute the conjugacy class of g
        current_class = set()
        for h in group:
            h_inv = invert_perm(h)
            # conjugate = h * g * h_inv
            conjugate = multiply_perms(h, multiply_perms(g, h_inv))
            current_class.add(conjugate)
            
        conjugacy_classes.append(frozenset(current_class))
        
        # Remove elements of this class from the remaining set
        group_elements.difference_update(current_class)
        
    return conjugacy_classes

def solve_problem():
    """
    Main function to solve the problem by calculating the number of
    conjugacy classes of A5.
    """
    print("Step 1: Determine the group G and its properties.")
    print("The orientation-preserving icosahedral group is isomorphic to A5, the alternating group on 5 elements.")
    print("The action is on C^3 and is in SL(3, C).")
    print("-" * 20)

    print("Step 2: Relate the cohomology to group theory via McKay Correspondence.")
    print("rank H^2_c(Y, Q) = b_2(Y) = (# non-trivial irreducible representations of A5).")
    print("This is equal to (# conjugacy classes of A5) - 1.")
    print("-" * 20)
    
    print("Step 3: Compute the number of conjugacy classes of A5.")
    
    # Generate S5
    s5 = list(itertools.permutations(range(5)))
    
    # Filter for A5 (even permutations)
    a5 = [p for p in s5 if get_permutation_sign(p) == 1]
    
    print(f"The order of A5 is |A5| = {len(a5)}.")
    
    # Calculate conjugacy classes
    classes = calculate_conjugacy_classes(a5)
    num_classes = len(classes)
    
    print(f"The elements of A5 are partitioned into {num_classes} conjugacy classes.")
    print("-" * 20)

    print("Step 4: Calculate the final answer.")
    num_nontrivial_irreps = num_classes - 1
    
    print(f"Number of irreducible representations = {num_classes}")
    print(f"Number of non-trivial irreducible representations = {num_classes} - 1 = {num_nontrivial_irreps}")
    print(f"Therefore, the rank of H^2_c(Y, Q) is {num_nontrivial_irreps}.")

if __name__ == '__main__':
    solve_problem()
<<<4>>>