import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def main():
    """
    Calculates the number of non-isomorphic functional graphs on 4 vertices.
    """
    print("This script calculates the number of conjugacy classes of endomorphisms on a set of size 4.")
    print("This is equivalent to counting non-isomorphic functional graphs on 4 vertices.")
    print("\nStep 1: Count connected functional graphs c(k) for k=1, 2, 3, 4.")

    # c(1): A single vertex must form a 1-cycle (a fixed point).
    c1 = 1
    print(f"c(1) = {c1}: A single vertex in a 1-cycle.")

    # c(2): On 2 vertices. Can be a 2-cycle or a 1-cycle with a tail.
    c2 = 1 + 1
    print(f"c(2) = {c2}: A 2-cycle, or a 1-cycle with a single vertex attached.")

    # c(3): On 3 vertices.
    # - A 3-cycle (1)
    # - A 2-cycle with a tail (1)
    # - A 1-cycle with a rooted tree on the other 2 vertices (number of rooted trees on 3 vertices is 2).
    c3 = 1 + 1 + 2
    print(f"c(3) = {c3}: A 3-cycle, a 2-cycle with a tail, or a 1-cycle with one of 2 possible rooted trees on 3 vertices.")

    # c(4): On 4 vertices.
    # - A 4-cycle (1)
    # - A 3-cycle with a tail (1)
    # - A 2-cycle with 2 other vertices forming trees rooted at the cycle (3 structures).
    # - A 1-cycle with a rooted tree on the other 3 vertices (number of rooted trees on 4 vertices is 4).
    c4 = 1 + 1 + 3 + 4
    print(f"c(4) = {c4}: A 4-cycle, a 3-cycle with a tail, one of 3 structures with a 2-cycle, or one of 4 rooted trees on 4 vertices.")

    print("\nStep 2: Sum contributions from each integer partition of 4.")

    # Partition 4: One component of size 4. Contribution = c(4)
    p4 = c4
    print(f"Partition 4: c(4) = {p4}")

    # Partition 3+1: One component of size 3, one of size 1. Contribution = c(3) * c(1)
    p31 = c3 * c1
    print(f"Partition 3+1: c(3) * c(1) = {c3} * {c1} = {p31}")

    # Partition 2+2: Two components of size 2. Choose 2 from c(2) types with replacement.
    # Contribution = C(c(2) + 2 - 1, 2) = C(3, 2) = 3
    p22 = combinations(c2 + 2 - 1, 2)
    print(f"Partition 2+2: C(c(2)+2-1, 2) = C({c2+1}, 2) = {p22}")

    # Partition 2+1+1: One of size 2, two of size 1. Contribution = c(2) * C(c(1) + 2 - 1, 2)
    p211 = c2 * combinations(c1 + 2 - 1, 2)
    print(f"Partition 2+1+1: c(2) * C(c(1)+2-1, 2) = {c2} * C({c1+1}, 2) = {p211}")

    # Partition 1+1+1+1: Four of size 1. Contribution = C(c(1) + 4 - 1, 4)
    p1111 = combinations(c1 + 4 - 1, 4)
    print(f"Partition 1+1+1+1: C(c(1)+4-1, 4) = C({c1+3}, 4) = {p1111}")

    total_count = p4 + p31 + p22 + p211 + p1111
    print("\nStep 3: The total number is the sum of these contributions.")
    print(f"Total = {p4} + {p31} + {p22} + {p211} + {p1111} = {total_count}")

if __name__ == "__main__":
    main()