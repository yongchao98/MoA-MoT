import math

def my_gcd(a, b):
    """A gcd function that handles the identity element 'inf'."""
    if a == 'inf':
        return b
    if b == 'inf':
        return a
    return math.gcd(a, b)

def main():
    """
    Calculates the cardinality of [F * G](.) for n=8.
    """
    n = 8
    # The set A is {1, ..., 8, inf}
    domain = list(range(1, n + 1)) + ['inf']
    domain_map = {val: i for i, val in enumerate(domain)}

    # Step 1: Generate the set of CMon homomorphisms from A to A.
    # We represent each function as a tuple of its values on the domain.
    homomorphisms = set()

    # Homomorphism 1: The constant map to infinity
    f_const_inf = tuple(['inf'] * len(domain))
    homomorphisms.add(f_const_inf)

    # Homomorphism 2: The identity map
    f_id = tuple(domain)
    homomorphisms.add(f_id)

    # Homomorphism family 3: f_d(x) = gcd(d, x) for x!=inf, f_d(inf)=inf
    f_d_family = []
    for d in range(1, n + 1):
        vals = [my_gcd(d, x) if x != 'inf' else 'inf' for x in domain]
        f_d_tuple = tuple(vals)
        f_d_family.append(f_d_tuple)
        homomorphisms.add(f_d_tuple)

    # Step 2: Calculate the cardinality of the tensor product.
    # This is the number of unique functions g(x) = f(gcd(b, x)).

    # We will build the set of unique functions and show how each family contributes.
    
    # Contribution from f_const_inf
    # g(x) = f_const_inf(gcd(b, x)) = 'inf' for all x. This gives one function.
    g_const_inf = tuple(['inf'] * len(domain))
    set1 = {g_const_inf}
    count1 = len(set1)
    
    # Contribution from f_id
    # g(x) = f_id(gcd(b, x)) = gcd(b, x). Let's call these h_b.
    # As b ranges over A, we get 9 distinct functions h_b, since h_b('inf') = b.
    set2 = set()
    for b in domain:
        vals = [my_gcd(b, x) for x in domain]
        set2.add(tuple(vals))
    count2 = len(set2)

    # Contribution from the f_d family
    # g(x) = f_d(gcd(b, x)).
    # If b != 'inf', g(x) = gcd(d, gcd(b, x)) = h_gcd(d,b)(x). These are already in set2.
    # If b == 'inf', g(x) = f_d(gcd('inf', x)) = f_d(x). This gives 8 new functions.
    set3 = set()
    for f_d in f_d_family:
        # This corresponds to b = 'inf'
        set3.add(f_d)
    
    # The sets are disjoint.
    # set1 vs set2: h_b('inf') = b != 'inf'. Disjoint.
    # set1 vs set3: f_d('inf') = 'inf', but f_d(d) = d != 'inf'. Disjoint.
    # set2 vs set3: h_b('inf') = b, f_d('inf') = 'inf'. Disjoint unless b='inf'.
    # h_inf(x) = x. f_d(x) = gcd(d,x). These are not equal for d in {1..8}.
    # So the three sets of functions are disjoint.
    
    count3 = len(set3)
    total_cardinality = count1 + count2 + count3

    print("The calculation of the cardinality is broken down as follows:")
    print(f"Number of unique functions generated from f_const_inf: {count1}")
    print(f"Number of unique functions generated from f_id: {count2}")
    print(f"Number of new unique functions generated from the f_d family: {count3}")
    print(f"The total cardinality is the sum of these disjoint sets.")
    print(f"Total cardinality = {count1} + {count2} + {count3} = {total_cardinality}")

if __name__ == "__main__":
    main()