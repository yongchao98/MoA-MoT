import itertools

def perm_compose(p2, p1):
    """Computes the composition of two permutations p2 . p1."""
    return tuple(p2[p1[i]] for i in range(len(p1)))

def main():
    """
    Finds all minimal grid diagrams for the left-hand trefoil knot.
    """
    n = 3
    # The canonical permutation for the left-hand trefoil in a normalized
    # (identity sigma) grid diagram.
    tau_L = (2, 0, 1)

    # Generate all possible permutations for sigma (the 'X' permutation).
    all_sigmas = list(itertools.permutations(range(n)))

    # The set of (sigma, tau) pairs that represent the left-hand trefoil.
    lh_trefoil_diagrams = []

    print("The minimal grid diagrams for the left-hand trefoil knot are represented by the following (sigma, tau) pairs:")
    print("---------------------------------------------------------------")
    
    for sigma in all_sigmas:
        # For a diagram (sigma, tau) to be a left-hand trefoil, its
        # normalized form must be (id, tau_L).
        # The normalization is tau_new = tau . sigma_inverse.
        # So, we need tau . sigma_inverse = tau_L.
        # This implies tau = tau_L . sigma.
        tau = perm_compose(tau_L, sigma)
        lh_trefoil_diagrams.append((sigma, tau))
        print(f"Diagram {len(lh_trefoil_diagrams)}: X-coords=({sigma[0]},{sigma[1]},{sigma[2]}), O-coords=({tau[0]},{tau[1]},{tau[2]})")

    print("---------------------------------------------------------------")
    count = len(lh_trefoil_diagrams)
    print(f"There are a total of {count} distinct minimal grid diagrams for the left-hand trefoil knot.")

if __name__ == "__main__":
    main()
