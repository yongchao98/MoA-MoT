import collections

def solve_quadratic_forms_z8():
    """
    Calculates the number of equivalence classes of quadratic forms in two
    variables over the ring R = Z/8Z.

    A quadratic form Q(x, y) = ax^2 + bxy + cy^2 is represented by a tuple (a, b, c).
    Two forms are equivalent if one can be obtained from the other by an
    invertible linear transformation of variables.

    This script iterates through all 512 possible forms, and for each unclassified
    form, it explores its entire equivalence class (orbit) using BFS.
    """
    N = 8

    # A list to keep track of visited forms, indexed by mapping (a,b,c) to
    # an integer key: a*N^2 + b*N + c.
    visited = [False] * (N * N * N)

    # First, generate all invertible 2x2 matrices over Z/8Z.
    # A matrix is invertible if its determinant is a unit in Z/8Z {1, 3, 5, 7}.
    invertible_matrices = []
    units = {1, 3, 5, 7}
    for p in range(N):
        for q_mat in range(N):
            for r in range(N):
                for s in range(N):
                    determinant = (p * s - q_mat * r) % N
                    if determinant in units:
                        invertible_matrices.append(((p, q_mat), (r, s)))

    num_classes = 0
    orbit_sizes = []

    # Iterate through all possible forms (a, b, c).
    for a_start in range(N):
        for b_start in range(N):
            for c_start in range(N):
                start_idx = a_start * N * N + b_start * N + c_start
                
                # If this form has not been visited, it's a new class.
                if not visited[start_idx]:
                    num_classes += 1
                    
                    # Begin BFS to find all forms in the class (the orbit).
                    q = collections.deque([(a_start, b_start, c_start)])
                    visited[start_idx] = True
                    current_orbit_size = 1
                    
                    while q:
                        a, b, c = q.popleft()

                        # Apply all transformations to the current form.
                        for p, r, s, q_mat in [(mat[0][0], mat[1][0], mat[1][1], mat[0][1]) for mat in invertible_matrices]:
                            # New coefficients (a', b', c') after transformation
                            a_new = (a*p*p + b*p*r + c*r*r) % N
                            b_new = (2*a*p*q_mat + b*(p*s + q_mat*r) + 2*c*r*s) % N
                            c_new = (a*q_mat*q_mat + b*q_mat*s + c*s*s) % N
                            
                            new_idx = a_new * N * N + b_new * N + c_new
                            
                            if not visited[new_idx]:
                                visited[new_idx] = True
                                q.append((a_new, b_new, c_new))
                                current_orbit_size += 1
                    
                    orbit_sizes.append(current_orbit_size)

    print(f"The total number of equivalence classes is: {num_classes}")
    orbit_sizes.sort()
    print(f"The sizes of the orbits are: {orbit_sizes}")
    equation_str = " + ".join(map(str, orbit_sizes))
    print(f"The sum of the orbit sizes is: {equation_str} = {sum(orbit_sizes)}")

if __name__ == '__main__':
    solve_quadratic_forms_z8()