import math

def find_solutions(max_n_to_check):
    """
    Finds and prints integer values of n up to a given limit for which the
    associate ring graph AG(Z_n) has an equal number of vertices and edges.

    This is based on the condition Sum[m_i * (m_i - 3)] = 0, where m_i are
    the sizes of the associate classes.
    """
    solution_values = []
    
    for n in range(3, max_n_to_check + 1):
        # The vertices are the non-zero elements {1, 2, ..., n-1}
        all_vertices = set(range(1, n))
        
        # The units are elements coprime to n
        units = {i for i in range(1, n) if math.gcd(i, n) == 1}
        
        class_sizes = []
        processed_vertices = set()
        
        # Partition the vertices into associate classes
        for i in range(1, n):
            if i in processed_vertices:
                continue
            
            # The associate class of i is the set {i*u mod n}
            current_class = {(i * u) % n for u in units}
            # The vertices are non-zero, so we discard 0 if it appears
            current_class.discard(0)

            # An element must be in its own class, so if the class is empty,
            # it means all associates resulted in 0. The class is just the element itself.
            if not current_class:
                current_class = {i}

            class_sizes.append(len(current_class))
            processed_vertices.update(current_class)
            
        # Check if the sum of class sizes matches the number of vertices.
        # This is a sanity check for the partition logic.
        if sum(class_sizes) != n - 1:
            # If the partitioning logic is complex, it might miscount.
            # A more robust method for partitioning:
            vertices_to_process = set(range(1,n))
            robust_sizes = []
            while vertices_to_process:
                v = vertices_to_process.pop()
                cls = {(v * u) % n for u in units}
                cls.discard(0)
                if not cls:
                    cls = {v}
                robust_sizes.append(len(cls))
                vertices_to_process.difference_update(cls)
            class_sizes = robust_sizes

        # Check the derived condition: Sum[m_i * (m_i - 3)] = 0
        condition_sum = sum(m * (m - 3) for m in class_sizes)
        
        if condition_sum == 0:
            solution_values.append(n)
            
    # Print the final result in the requested format.
    # The phrase "output each number in the final equation" is interpreted as
    # showing the calculation for the found solution(s).
    if not solution_values:
        print(f"No solutions found for n up to {max_n_to_check}.")
    else:
        for n_sol in solution_values:
            # For showing the equation, we re-calculate the classes for the solution
            units = {i for i in range(1, n_sol) if math.gcd(i, n_sol) == 1}
            vertices_to_process = set(range(1,n_sol))
            final_class_sizes = []
            while vertices_to_process:
                v = vertices_to_process.pop()
                cls = {(v * u) % n_sol for u in units}
                cls.discard(0)
                if not cls: cls = {v}
                final_class_sizes.append(len(cls))
                vertices_to_process.difference_update(cls)

            print(f"Solution found: n = {n_sol}")
            print(f"  The associate class sizes are: {sorted(final_class_sizes)}")
            eq_parts = [f"{m}({m}-3)" for m in sorted(final_class_sizes)]
            eq_sum = sum(m * (m - 3) for m in final_class_sizes)
            print(f"  The final equation is: {' + '.join(eq_parts)} = {eq_sum}")

        sequence = ", ".join(map(str, solution_values))
        print(f"\nThus, the set of values is n \u2208 {{ {sequence} }}")

# Run the search for n up to a reasonable limit, e.g., 1000.
find_solutions(1000)