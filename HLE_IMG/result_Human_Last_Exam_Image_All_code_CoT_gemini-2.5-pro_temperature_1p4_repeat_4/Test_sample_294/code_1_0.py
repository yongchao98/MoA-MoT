def find_minimal_N():
    """
    Searches for the smallest N for which S(N) is non-empty.
    S(N) contains pairs (n1, n2) with n1,n2 <= N that can generate an infinite group.
    
    The condition for an infinite group is n1>=2, n2>=2, and 1/n1 + 1/n2 <= 1.
    """
    N = 1
    while True:
        print(f"Checking for N = {N}...")
        found_pair = None
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # An infinite group requires n1 >= 2 and n2 >= 2.
                if n1 >= 2 and n2 >= 2:
                    # Check if the pair corresponds to a Euclidean or Hyperbolic case.
                    if (1/n1 + 1/n2) <= 1:
                        found_pair = (n1, n2)
                        break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            print(f"For N = {N}, found the pair (n1, n2) = ({n1}, {n2}).")
            print("This pair satisfies n1 <= N and n2 <= N.")
            print("Checking the condition for this pair to generate an infinite group:")
            result = 1/n1 + 1/n2
            print(f"1 / {n1} + 1 / {n2} = {result}")
            if result == 1:
                print(f"The result is equal to 1 (Euclidean case), so the group can be infinite.")
            else:
                 print(f"The result is less than 1 (Hyperbolic case), so the group can be infinite.")
            print(f"\nSince a valid pair was found for N = {N}, and not for smaller N, this is the minimum value.")
            return N
        else:
            print(f"No pair (n1, n2) satisfying the conditions found for N = {N}. S({N}) is empty.")
            N += 1

min_N = find_minimal_N()
print(f"\nThe solution is {min_N}.")