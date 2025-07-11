import math

def factorial(n):
    """Computes factorial of n."""
    if n == 0:
        return 1
    return math.factorial(n)

def calculate_S():
    """
    Calculates the total number of ways to distribute the items (S).
    S = 25! / (5!)^5
    """
    return factorial(25) // (factorial(5)**5)

def W(block_matrix):
    """
    Calculates the number of arrangements for a given block matrix.
    W(C) = (5!)^n / product(Cij!) where n is the size of the block.
    """
    n = len(block_matrix)
    numerator = factorial(5)**n
    denominator = 1
    for i in range(n):
        for j in range(n):
            denominator *= factorial(block_matrix[i][j])
    return numerator // denominator

def calculate_F_fixed():
    """
    Calculates the number of favorable distributions for a fixed specialization.
    This is done by summing contributions from all partitions of the 5 individuals.
    """
    # Calculate arrangements for irreducible blocks I(k)
    I = {}
    # I(1): block is [[5]]
    I[1] = 1
    
    # I(2): blocks are [[4,1],[1,4]] and [[3,2],[2,3]]
    I[2] = W([[4, 1], [1, 4]]) + W([[3, 2], [2, 3]])

    # I(3): asymmetric block (3 ways to choose special index) and symmetric block
    b3_asym = W([[3, 1, 1], [1, 4, 0], [1, 0, 4]])
    b3_sym = W([[3, 1, 1], [1, 3, 1], [1, 1, 3]])
    I[3] = 3 * b3_asym + b3_sym

    # I(4): asymmetric block (4 ways to choose special index)
    b4_asym = W([[2, 1, 1, 1], [1, 4, 0, 0], [1, 0, 4, 0], [1, 0, 0, 4]])
    I[4] = 4 * b4_asym

    # I(5): No irreducible 5x5 blocks exist for this problem
    I[5] = 0

    # Sum contributions from all partitions of the set {1,2,3,4,5}
    f_fixed = 0
    # Partition {1,1,1,1,1}
    f_fixed += I[1]**5
    # Partition {2,1,1,1}
    f_fixed += (factorial(5) // (factorial(2) * factorial(3))) * I[2] * (I[1]**3)
    # Partition {2,2,1}
    f_fixed += (factorial(5) // (factorial(2) * factorial(2) * factorial(1) * factorial(2))) * (I[2]**2) * I[1]
    # Partition {3,1,1}
    f_fixed += (factorial(5) // (factorial(3) * factorial(2))) * I[3] * (I[1]**2)
    # Partition {3,2}
    f_fixed += (factorial(5) // (factorial(3) * factorial(2))) * I[3] * I[2]
    # Partition {4,1}
    f_fixed += (factorial(5) // (factorial(4) * factorial(1))) * I[4] * I[1]
    # Partition {5}
    f_fixed += I[5]
    
    return f_fixed

def main():
    """
    Main function to calculate S, F, and the final probability P.
    """
    S = calculate_S()
    F_fixed = calculate_F_fixed()
    F = factorial(5) * F_fixed
    P = F / S

    print(f"Total number of distributions (S): {S}")
    print(f"Number of favorable distributions (F): {F}")
    print(f"The probability P is the ratio F/S.")
    print(f"P = F/S = {F}/{S} = {P}")

if __name__ == "__main__":
    main()