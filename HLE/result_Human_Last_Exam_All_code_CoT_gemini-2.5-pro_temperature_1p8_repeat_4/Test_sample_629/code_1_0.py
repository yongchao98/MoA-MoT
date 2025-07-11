import itertools

def perm_inverse(p):
    """Computes the inverse of a permutation."""
    inv = [0] * len(p)
    for i, val in enumerate(p):
        inv[val] = i
    return tuple(inv)

def perm_compose(p1, p2):
    """Computes the composition of two permutations, p1 o p2."""
    return tuple(p1[p2[i]] for i in range(len(p1)))

def apply_translation(diagram, a, b, n=3):
    """Applies a translation by (a, b) to a diagram."""
    sigma, pi = diagram
    # Column shift by a
    sigma_col_shifted = tuple(sigma[(i - a) % n] for i in range(n))
    pi_col_shifted = tuple(pi[(i - a) % n] for i in range(n))
    # Row shift by b
    new_sigma = tuple((val + b) % n for val in sigma_col_shifted)
    new_pi = tuple((val + b) % n for val in pi_col_shifted)
    return (new_sigma, new_pi)

def main():
    """Main function to perform the calculation and print the results."""
    n = 3
    perms = list(itertools.permutations(range(n)))

    # The left-hand trefoil corresponds to one of the two 3-cycles.
    LEFT_TREFOIL_PERM = (1, 2, 0) # This represents c(i) = i+1 mod 3

    # 1. Generate all grid diagrams for the left-hand trefoil.
    left_trefoil_diagrams = []
    for sigma in perms:
        for pi in perms:
            pi_inv = perm_inverse(pi)
            composed_perm = perm_compose(pi_inv, sigma)
            if composed_perm == LEFT_TREFOIL_PERM:
                left_trefoil_diagrams.append((sigma, pi))
    
    print(f"Found {len(left_trefoil_diagrams)} minimal grid diagrams for the left-hand trefoil.\n")

    # 2. Find the number of orbits under the translation group action.
    unclassified_diagrams = set(left_trefoil_diagrams)
    orbit_count = 0
    orbits = []
    
    while unclassified_diagrams:
        orbit_count += 1
        start_diagram = unclassified_diagrams.pop()
        current_orbit = {start_diagram}
        queue = [start_diagram]
        head = 0
        while head < len(queue):
            current_diagram = queue[head]
            head += 1
            for a in range(n):
                for b in range(n):
                    translated_diagram = apply_translation(current_diagram, a, b, n)
                    if translated_diagram not in current_orbit:
                        current_orbit.add(translated_diagram)
                        queue.append(translated_diagram)
        unclassified_diagrams -= current_orbit
        orbits.append(sorted(list(current_orbit)))
    
    print(f"These {len(left_trefoil_diagrams)} diagrams fall into {orbit_count} orbits under translation:")
    for i, orbit in enumerate(orbits):
        print(f"\nOrbit {i+1}:")
        for d_sigma, d_pi in orbit:
            print(f"  σ={d_sigma}, π={d_pi}")
    
    print("\n-------------------------")
    print("Final Calculation Summary")
    print("-------------------------")
    
    # 3. Explain the result using Burnside's Lemma as a final check.
    #    Number of orbits = (1 / |Group Size|) * Σ |Diagrams fixed by g|
    
    fixed_point_counts = {}
    for a in range(n):
        for b in range(n):
            g = (a, b)
            count = 0
            for diagram in left_trefoil_diagrams:
                if apply_translation(diagram, a, b) == diagram:
                    count += 1
            fixed_point_counts[g] = count
    
    sum_of_fixed_points = sum(fixed_point_counts.values())

    print("To formally count the number of distinct diagrams under translation, we use Burnside's Lemma.")
    print("Number of orbits = (Sum of diagrams fixed by each translation) / (Total number of translations)")
    print("\nSum of fixed diagrams:")
    fixed_values = list(fixed_point_counts.values())
    print(" + ".join(map(str, fixed_values)) + f" = {sum_of_fixed_points}")
    print("\nTotal number of translations on a 3x3 grid:")
    print("3 (column shifts) * 3 (row shifts) = 9")

    final_answer = sum_of_fixed_points // (n * n)
    print("\nFinal Equation:")
    print(f"{sum_of_fixed_points} / 9 = {final_answer}")
    
    print("\nAs explained, rotating a diagram for a left-hand trefoil yields a diagram for a right-hand trefoil.")
    print("Therefore, considering rotations does not reduce the number of distinct diagrams further.")
    print(f"\nThe total number of unique grid diagrams is {final_answer}.")

if __name__ == '__main__':
    main()