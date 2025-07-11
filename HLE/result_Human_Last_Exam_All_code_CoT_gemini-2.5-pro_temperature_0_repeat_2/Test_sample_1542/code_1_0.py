import numpy as np

def get_gl2_z8():
    """
    Generates the group GL_2(Z/8Z) of invertible 2x2 matrices over Z/8Z.
    A matrix is invertible if its determinant is a unit in Z/8Z (i.e., odd).
    """
    G = []
    for p in range(8):
        for q in range(8):
            for r in range(8):
                for s in range(8):
                    det = (p * s - q * r) % 8
                    if det in [1, 3, 5, 7]:
                        G.append(np.array([[p, q], [r, s]], dtype=int))
    return G

def apply_transform(T, Q):
    """
    Applies a linear transform T to a quadratic form Q.
    If Q is ax^2 + bxy + cy^2 and the transform is (x,y) -> (px+qy, rx+sy),
    the new form Q' = a'x^2 + b'xy + c'y^2 has coefficients:
    a' = ap^2 + bpr + cr^2
    b' = 2apq + b(ps+qr) + 2crs
    c' = aq^2 + bqs + cs^2
    """
    a, b, c = Q
    p, q = T[0, 0], T[0, 1]
    r, s = T[1, 0], T[1, 1]
    
    a_new = (a * p * p + b * p * r + c * r * r) % 8
    b_new = (2 * a * p * q + b * (p * s + q * r) + 2 * c * r * s) % 8
    c_new = (a * q * q + b * q * s + c * s * s) % 8
    
    return (a_new, b_new, c_new)

def count_classes_in_set(form_set, G):
    """
    Counts the number of equivalence classes (orbits) in a given set of forms.
    """
    num_classes = 0
    forms_to_process = form_set.copy()
    
    while forms_to_process:
        num_classes += 1
        # Pick a representative form from the set
        q_rep = forms_to_process.pop()
        
        # Compute its orbit
        orbit = set()
        queue = [q_rep]
        visited = {q_rep}
        
        head = 0
        while head < len(queue):
            current_q = queue[head]
            head += 1
            orbit.add(current_q)
            
            for T in G:
                new_q = apply_transform(T, current_q)
                if new_q not in visited:
                    visited.add(new_q)
                    queue.append(new_q)
        
        # Remove the entire orbit from the set of forms to process
        forms_to_process -= orbit
        
    return num_classes

def solve():
    """
    Main function to classify quadratic forms and print the result.
    """
    # Partition all 512 forms into three disjoint sets based on their type.
    forms_type_I_b_odd = set()
    forms_type_I_b_even = set()
    forms_type_II = set()

    for a in range(8):
        for b in range(8):
            for c in range(8):
                q = (a, b, c)
                # Type II: a, b, c are all even
                if a % 2 == 0 and b % 2 == 0 and c % 2 == 0:
                    forms_type_II.add(q)
                # Type I, b odd
                elif b % 2 != 0:
                    forms_type_I_b_odd.add(q)
                # Type I, b even (and not all a,b,c are even)
                else:
                    forms_type_I_b_even.add(q)

    # Generate the transformation group
    G = get_gl2_z8()

    # Count classes in each partition
    count_I_b_even = count_classes_in_set(forms_type_I_b_even, G)
    count_I_b_odd = count_classes_in_set(forms_type_I_b_odd, G)
    count_II = count_classes_in_set(forms_type_II, G)
    
    total = count_I_b_even + count_I_b_odd + count_II
    
    print("The number of equivalence classes can be found by classifying the forms into types.")
    print(f"Number of classes of Type I (b even): {count_I_b_even}")
    print(f"Number of classes of Type I (b odd): {count_I_b_odd}")
    print(f"Number of classes of Type II (a,b,c all even): {count_II}")
    print("\nThe total number of equivalence classes is the sum of these counts:")
    print(f"{count_I_b_even} + {count_I_b_odd} + {count_II} = {total}")

solve()
<<<42>>>