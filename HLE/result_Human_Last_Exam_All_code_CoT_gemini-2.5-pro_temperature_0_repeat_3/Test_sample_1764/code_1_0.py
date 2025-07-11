import itertools

def demonstrate_minimum_embeddings():
    """
    This script demonstrates that the smallest possible number of isometric embeddings
    from a finite ultrametric space X into a Banach space B is 0.
    """

    # 1. Define a simple finite ultrametric space X.
    # Let X = {a, b} with d(a, b) = 1.
    # This is an ultrametric space because the strong triangle inequality holds for all triplets.
    X_points = ['a', 'b']
    
    def distance_in_X(p1, p2):
        if p1 == p2:
            return 0
        else:
            return 1

    # 2. Define the simplest Banach space B.
    # Let B = {0}, the trivial Banach space containing only the zero vector.
    # We represent the vector 0 as the number 0 for simplicity.
    B_vectors = [0]
    
    def norm_in_B(v):
        # The norm of the zero vector is 0.
        return abs(v)

    # 3. Search for all possible isometric embeddings from X to B.
    
    # An embedding is a function f: X -> B. Let's generate all possible functions.
    # In this case, there is only one possible function: f(a)=0, f(b)=0.
    
    num_embeddings = 0
    
    # Iterate through all possible mappings from X_points to B_vectors
    for p_images in itertools.product(B_vectors, repeat=len(X_points)):
        f = dict(zip(X_points, p_images))
        
        is_isometry = True
        # Check the isometry condition: ||f(p1) - f(p2)|| = d(p1, p2)
        for p1 in X_points:
            for p2 in X_points:
                dist_X = distance_in_X(p1, p2)
                
                # In B={0}, vector subtraction is just 0 - 0 = 0.
                dist_B = norm_in_B(f[p1] - f[p2])
                
                if dist_X != dist_B:
                    is_isometry = False
                    # We found a pair of points for which the distance is not preserved.
                    # Let's print the details of this contradiction.
                    print("Checking mapping f where f(a)={} and f(b)={}...".format(f['a'], f['b']))
                    print("For points '{}' and '{}':".format(p1, p2))
                    print("  - Distance in X: d({}, {}) = {}".format(p1, p2, dist_X))
                    print("  - Distance in B: ||f({}) - f({})|| = {}".format(p1, p2, dist_B))
                    print("The equation d(p1, p2) = ||f(p1) - f(p2)|| becomes:")
                    print("  {} = {}".format(dist_X, dist_B))
                    print("This is false, so this mapping is not an isometry.")
                    break
            if not is_isometry:
                break
        
        if is_isometry:
            num_embeddings += 1

    print("\n--- Conclusion ---")
    print(f"The total number of isometric embeddings found is: {num_embeddings}")
    print("This example shows that it is possible for the number of embeddings to be 0.")
    print("Since the number of embeddings cannot be negative, the smallest possible number is 0.")

if __name__ == '__main__':
    demonstrate_minimum_embeddings()
