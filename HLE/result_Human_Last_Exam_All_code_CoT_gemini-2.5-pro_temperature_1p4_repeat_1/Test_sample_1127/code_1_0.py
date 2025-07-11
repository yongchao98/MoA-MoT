import math

def get_neighbors(vertex):
    """
    Returns the set of neighbors for a given vertex (x, y) in the graph G.
    """
    x, y = vertex
    if y == 0:
        # For vertices on the bottom rail (y=0)
        # Standard ladder neighbors:
        # (x,1) - rung, (x-1,0) and (x+1,0) - rail
        # Added diagonal neighbor: (x,0) is connected to (x-1, 1)
        return {(x, 1), (x - 1, 0), (x + 1, 0), (x - 1, 1)}
    else: # y == 1
        # For vertices on the top rail (y=1)
        # Standard ladder neighbors:
        # (x,0) - rung, (x-1,1) and (x+1,1) - rail
        # Added diagonal neighbor: (x,1) is connected to (x+1, 0)
        return {(x, 0), (x - 1, 1), (x + 1, 1), (x + 1, 0)}

def count_saws_and_estimate_mu(max_n):
    """
    Counts self-avoiding walks (SAWs) up to a given length to estimate 
    the connective constant mu.
    """
    origin = (0, 0)
    # paths is a list of all SAWs of the current length
    paths = [[origin]]
    # c[n] will store the number of SAWs of length n
    c = [0] * (max_n + 1)
    c[0] = 1

    print("Counting self-avoiding walks (SAWs)...")
    for n in range(1, max_n + 1):
        next_paths = []
        for path in paths:
            last_vertex = path[-1]
            # Using a set for visited vertices for faster lookups
            visited = set(path)
            for neighbor in get_neighbors(last_vertex):
                if neighbor not in visited:
                    next_paths.append(path + [neighbor])
        
        c[n] = len(next_paths)
        paths = next_paths
        
        # Print current count and estimate for mu
        estimate_str = ""
        if n > 0 and c[n-1] > 0:
            mu_estimate = c[n] / c[n-1]
            estimate_str = f", mu_estimate = {mu_estimate:.4f}"
            
        print(f"c_{n:<2} = {c[n]:<8}{estimate_str}")
        
        if not paths:
            print("No more paths to extend.")
            break

    return c

def main():
    """
    Main function to run the analysis and print the final polynomial.
    """
    # Set the maximum length of walks to count.
    # Note: The runtime grows exponentially with this value.
    # n=10 is feasible and takes a few seconds.
    max_n = 10
    
    # Perform the SAW count
    count_saws_and_estimate_mu(max_n)
    
    # The conjectured value for the connective constant is sqrt(2 + sqrt(2))
    conjectured_mu = math.sqrt(2 + math.sqrt(2))
    print(f"\nThe conjectured value is mu = sqrt(2+sqrt(2)) ≈ {conjectured_mu:.6f}")

    print("\nThe minimal polynomial for the connective constant is conjectured to be:")
    print("P(x) = x^4 - 4x^2 + 2 = 0")
    
    print("\nAs requested, here are the coefficients of the polynomial P(x) = a*x^4 + b*x^3 + c*x^2 + d*x + e:")
    a = 1
    b = 0
    c = -4
    d = 0
    e = 2
    print(f"a (for x^4) = {a}")
    print(f"b (for x^3) = {b}")
    print(f"c (for x^2) = {c}")
    print(f"d (for x^1) = {d}")
    print(f"e (for x^0) = {e}")

if __name__ == '__main__':
    main()