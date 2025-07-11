import math

def get_phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

def get_associate_class_sizes(n):
    """
    Computes the sizes of all associate classes for Z_n.
    An associate class of element a is the set {a*u mod n for u in units}.
    """
    if n <= 1:
        return []

    # Find units
    units = [i for i in range(1, n) if math.gcd(i, n) == 1]
    
    non_zero_elements = set(range(1, n))
    classified_elements = set()
    class_sizes = []

    for i in sorted(list(non_zero_elements)):
        if i not in classified_elements:
            current_class = set()
            for u in units:
                current_class.add((i * u) % n)
            
            # Remove 0 if it was generated (shouldn't happen for i in non-zero-divisors)
            if 0 in current_class:
                current_class.remove(0)

            class_sizes.append(len(current_class))
            classified_elements.update(current_class)
            
    return class_sizes

def find_n_for_ring_graph():
    """
    Finds values of n for which AG(Z_n) is a "ring graph",
    interpreted as a graph where all components are paths (K1 or K2).
    This means all associate class sizes must be 1 or 2.
    """
    solutions = []
    # We test n up to a reasonable limit to find the sequence.
    # The theoretical argument shows we only need to check where phi(n)<=2.
    # Those are n=2,3,4,6. We check a bit beyond to be sure.
    for n in range(2, 50):
        # A necessary condition is that the class of units has size <= 2
        if get_phi(n) > 2:
            continue
        
        # Now verify all other classes also meet the condition
        class_sizes = get_associate_class_sizes(n)
        is_ring_graph = all(size <= 2 for size in class_sizes)
        
        if is_ring_graph:
            solutions.append(n)
            
    # As per convention, n=2 (which yields a single vertex graph) is often excluded
    final_solutions = [n for n in solutions if n > 2]
    return final_solutions

# Calculate and print the result
result_values = find_n_for_ring_graph()
print(f"The values of n for which AG(Z_n) is a ring graph are n in {{{', '.join(map(str, result_values))}}}.")
