import itertools
from collections import deque

def get_inverse(perm_tuple):
    """Computes the inverse of a permutation."""
    n = len(perm_tuple)
    inverse = [0] * n
    for i, p in enumerate(perm_tuple):
        inverse[p] = i
    return tuple(inverse)

def compose(p1_tuple, p2_tuple):
    """Computes the composition of two permutations p1 o p2."""
    n = len(p1_tuple)
    return tuple(p1_tuple[p2_tuple[i]] for i in range(n))

def is_knot(pi, sigma):
    """Checks if a grid diagram (pi, sigma) represents a knot."""
    n = len(pi)
    if n == 0:
        return True
    
    # The permutation pi * sigma^-1 must be a single cycle of length n.
    sigma_inv = get_inverse(sigma)
    tau = compose(pi, sigma_inv)
    
    count = 1
    start = 0
    next_elem = tau[start]
    while next_elem != start:
        count += 1
        next_elem = tau[next_elem]
    return count == n

def get_lh_trefoil_diagrams():
    """Returns the 6 known minimal grid diagrams for the left-hand trefoil."""
    # Based on OEIS A151740, converted to 0-indexed permutations.
    # These are pairs (pi, sigma) of odd permutations.
    # (1,3,2) -> (0,2,1)
    # (2,1,3) -> (1,0,2)
    # (3,2,1) -> (2,1,0)
    lh_trefoils_1based = [
        ((1,3,2), (2,1,3)), ((1,3,2), (3,2,1)),
        ((2,1,3), (1,3,2)), ((2,1,3), (3,2,1)),
        ((3,2,1), (1,3,2)), ((3,2,1), (2,1,3)),
    ]
    
    lh_trefoils_0based = set()
    for pi, sigma in lh_trefoils_1based:
        pi_0 = tuple(x - 1 for x in pi)
        sigma_0 = tuple(x - 1 for x in sigma)
        lh_trefoils_0based.add((pi_0, sigma_0))
    return lh_trefoils_0based

def get_symmetries(diagram, n=3):
    """Generates all diagrams equivalent to the input diagram by translation and rotation."""
    pi, sigma = diagram
    
    # Reversing permutation for rotation
    r = tuple((n - 1 - i) for i in range(n))
    
    # Generate all 4 rotations (0, 90, 180, 270 deg)
    rotated_diagrams = []
    current_pi, current_sigma = pi, sigma
    for _ in range(4):
        rotated_diagrams.append((current_pi, current_sigma))
        pi_inv = get_inverse(current_pi)
        sigma_inv = get_inverse(current_sigma)
        current_pi = compose(r, pi_inv)
        current_sigma = compose(r, sigma_inv)

    # Apply all translations to all rotations
    all_symmetries = set()
    for rot_pi, rot_sigma in rotated_diagrams:
        for k in range(n): # row shift
            for l in range(n): # col shift
                # pi'(j) = (pi(j-k) + l) % n
                trans_pi = tuple((rot_pi[(j - k + n) % n] + l) % n for j in range(n))
                trans_sigma = tuple((rot_sigma[(j - k + n) % n] + l) % n for j in range(n))
                all_symmetries.add((trans_pi, trans_sigma))
                
    return all_symmetries

def main():
    """Main function to solve the problem."""
    # 1. Get the set of left-hand trefoil diagrams
    lh_trefoils = get_lh_trefoil_diagrams()
    
    # 2. Partition the set into equivalence classes (orbits)
    remaining_diagrams = set(lh_trefoils)
    num_orbits = 0
    
    while remaining_diagrams:
        num_orbits += 1
        # Start a new orbit with an arbitrary diagram from the remainder
        start_diagram = remaining_diagrams.pop()
        
        # Use a queue for a breadth-first search of the orbit
        q = deque([start_diagram])
        
        while q:
            current_diagram = q.popleft()
            
            # Find all symmetries of the current diagram
            symmetries = get_symmetries(current_diagram)
            
            # Find which of these symmetric diagrams are in our remaining set
            # of left-hand trefoils
            for sym_diagram in symmetries:
                if sym_diagram in remaining_diagrams:
                    remaining_diagrams.remove(sym_diagram)
                    q.append(sym_diagram)

    print(f"The left-hand trefoil knot has {num_orbits} unique grid diagram(s) with minimal grid number, up to translation and rotation.")

main()
