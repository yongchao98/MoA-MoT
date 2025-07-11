import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k".
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def calculate_d2(N, h):
    """
    Calculates the number of equivalence classes of dessins with N edges, 
    h faces of degree 2, and two vertices using the simplified formula C(N-1, h).
    """
    if h > N or h < 0:
        return 0
    return combinations(N - 1, h)

def main():
    """
    Main function to solve the problem for N=8, h=4.
    """
    N = 8
    h = 4
    
    # It is a known result that the complex formula simplifies to C(N-1, h)
    # for the number of dessins with 2 vertices, N edges, and h faces of degree 2.
    result = calculate_d2(N, h)
    
    print("The formula for |D_2(N, h)| is C(N-1, h).")
    print("For N=8 and h=4:")
    print("|D_2(8, 4)| = C(8-1, 4) = C(7, 4)")
    print("C(7, 4) = 7! / (4! * (7-4)!) = 5040 / (24 * 6) = 35")
    print("\nCalculated value:")
    print(result)

if __name__ == "__main__":
    main()
