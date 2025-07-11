def analyze_complexity():
    """
    This function prints a step-by-step analysis of the given sorting
    algorithm's worst-case time complexity.
    """
    print("### Algorithm Time Complexity Analysis ###")

    print("\n1. Recurrence Relation Formulation:")
    print("Let T(n) be the time complexity for sorting n elements.")
    print("In the worst case, the algorithm performs two recursive calls on a subproblem of size n-1, plus a constant amount of work (c) for the comparison and swap.")
    print("This gives the recurrence relation: T(n) = 2 * T(n-1) + c")
    print("The base case is T(1) = c1 (constant time).")

    print("\n2. Solving the Recurrence Relation by Expansion:")
    print("T(n) = 2*T(n-1) + c")
    print("     = 2*(2*T(n-2) + c) + c = 4*T(n-2) + 2*c + c")
    print("     = 4*(2*T(n-3) + c) + 3*c = 8*T(n-3) + 4*c + 3*c")
    print("The general form after k steps is: T(n) = 2^k * T(n-k) + (2^k - 1)*c")
    
    print("\n3. Reaching the Base Case:")
    print("We set k = n-1 to get T(1) on the right side.")
    print("T(n) = 2^(n-1) * T(1) + (2^(n-1) - 1)*c")
    print("Substituting T(1) = c1:")
    print("T(n) = c1 * 2^(n-1) + c * 2^(n-1) - c")
    print("This can be written as: T(n) = (c1 + c) * 2^(n-1) - c")

    print("\n4. Final Complexity:")
    print("The dominant term is 2^(n-1). Therefore, the time complexity is O(2^n).")

if __name__ == '__main__':
    analyze_complexity()