def solve_modulo_permutation(initial_x, a):
    """
    Finds the optimal permutation of list 'a' to minimize the change in 'x'
    after sequential modulo operations.

    Prints the optimal permutation and the step-by-step calculation.
    """
    if not a:
        print(f"Input list 'a' is empty. No operations to perform.")
        print(f"Initial x = {initial_x}")
        print(f"Final x = {initial_x}")
        return

    # Find the minimum element in the list.
    try:
        a_min = min(a)
    except ValueError:
        print("Input list 'a' is empty.")
        return

    # If x is already smaller than the smallest element, its value will never change.
    # The absolute difference is 0, which is the minimum possible.
    if initial_x < a_min:
        print(f"Initial x = {initial_x}")
        # Any permutation is optimal. We can just use the original list 'a'.
        print(f"The best permutation is any permutation, for example: {a}")
        print(f"Since {initial_x} is smaller than the minimum element {a_min}, the value of x will not change.")
        
        equation_str = str(initial_x)
        for val in a:
            equation_str = f"({equation_str} mod {val})"
        print(f"The final equation is: {equation_str} = {initial_x}")
        print(f"Final x = {initial_x}")
        return

    # Create a list of elements other than a_min.
    a_prime = a.copy()
    a_prime.remove(a_min)

    # Candidate 1: permutation starts with a_min.
    max_x = initial_x % a_min
    best_b = None

    # Find the best candidate by iterating through all other elements `b` as
    # the first effective operator.
    for b in a_prime:
        candidate_x = (initial_x % b) % a_min
        if candidate_x > max_x:
            max_x = candidate_x
            best_b = b

    # Construct the optimal permutation based on our finding.
    p = []
    a_rem = a.copy()
    if best_b is None:
        # The best strategy is to apply a_min first.
        p.append(a_min)
        a_rem.remove(a_min)
        p.extend(a_rem)
    else:
        # The best strategy involves applying best_b, then a_min.
        p.append(best_b)
        a_rem.remove(best_b)
        p.append(a_min)
        a_rem.remove(a_min)
        p.extend(a_rem)

    # Print the results, including the step-by-step calculation.
    print(f"Initial x = {initial_x}")
    print(f"The best permutation found is: {p}")
    
    current_x = initial_x
    equation_str = str(initial_x)
    
    print("\nCalculation steps:")
    for val in p:
        next_x = current_x % val
        print(f"{current_x} mod {val} = {next_x}")
        current_x = next_x
        equation_str = f"({equation_str} mod {val})"

    final_x = current_x
    
    print(f"\nThe final equation is: {equation_str} = {final_x}")
    print(f"Final x = {final_x}")


# Example usage:
# x = 100, a = [13, 7] -> perm=[13,7] or [7,13], final_x=2
# x = 10, a = [6, 4] -> perm=[4,6], final_x=2
# x = 17, a = [12, 7] -> perm=[12,7], final_x=5

# Provide an example for execution
# You can change these values to test with other inputs
start_x = 17
numbers = [12, 7, 15]

solve_modulo_permutation(start_x, numbers)