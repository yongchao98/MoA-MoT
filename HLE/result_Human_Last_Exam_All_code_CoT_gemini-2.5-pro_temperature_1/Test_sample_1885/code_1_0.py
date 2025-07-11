# This script illustrates the concept of diagonalization, which is a key
# technique in constructing counterexamples in set theory. The user's
# question concerns objects (omega_1, omega_2) that cannot be represented
# or computed in Python, so this code is a conceptual illustration using
# finite objects and natural numbers.

def main():
    """
    Main function to run the demonstration.
    """
    # Let's imagine our functions map a finite domain to natural numbers.
    # In the actual problem, the domain is the uncountable set omega_1.
    domain_size = 10
    domain = range(domain_size)

    # In the actual proof (e.g., under CH), one must defeat all possible
    # bounding functions. Here, for demonstration, we create a small list
    # of "candidate" bounding functions that we want to diagonalize against.
    def g0(n):
        return n * n

    def g1(n):
        return 2**n + 5

    def g2(n):
        return 100

    candidate_bounds = [g0, g1, g2]
    
    print("Demonstrating diagonalization to construct an 'unbounded' sequence of functions.")
    print("-----------------------------------------------------------------------------")

    # We will construct a sequence of functions f_0, f_1, f_2, ...
    # Each function f_i will be constructed to "defeat" the i-th
    # candidate bounding function g_i. The sequence will also be strictly
    # increasing, analogous to the <* order in the problem.

    # We start with the zero function as our base case.
    f_previous = lambda n: 0
    
    # Iterate through the candidate bounds and construct a new function for each.
    for i, g in enumerate(candidate_bounds):
        # We define a new function f_new that "beats" the candidate bound g
        # and is greater than the previously constructed function f_previous.
        # The diagonalization happens at the coordinate 'i % domain_size'.
        
        # We use a closure to correctly capture the variables for f_new.
        def make_f_new(current_i, candidate_g, prev_f):
            def f_new(n):
                # Get the value from the previous function in the sequence
                prev_val = prev_f(n)
                
                # Diagonalization step: at a specific coordinate, we ensure f_new's
                # value is larger than the candidate bound g.
                diag_point = current_i % domain_size
                if n == diag_point:
                    g_val = candidate_g(n)
                    return max(prev_val, g_val) + 1
                else:
                    # For all other coordinates, just ensure the function value increases.
                    return prev_val + 1
            return f_new

        f_new = make_f_new(i, g, f_previous)
        
        # Update f_previous for the next iteration
        f_previous = f_new

        # --- Print the results of this construction step ---
        print(f"Step {i}:")
        print(f"  Constructed f_{i} to defeat the candidate bound g_{i} (named '{g.__name__}').")
        
        # Verify that f_i is indeed not bounded by g_i by finding a point
        # where f_i(n) >= g_i(n). Our construction guarantees this at 'i % domain_size'.
        unbounded_at = i % domain_size
        f_val = f_new(unbounded_at)
        g_val = g(unbounded_at)

        print(f"  f_{i} is not bounded by g_{i}. Check at n = {unbounded_at}:")
        print(f"    f_{i}({unbounded_at}) = {f_val}")
        print(f"    g_{i}({unbounded_at}) = {g_val}")
        
        # Output the comparison as an "equation" as requested.
        print(f"    Equation: {f_val} >= {g_val}")
        print("-" * 20)

    print("\nFinal constructed sequence has been created.")
    print("Each function f_i in the sequence is not bounded by the corresponding candidate g_i.")
    print("A transfinite version of this diagonalization process is used to build a counterexample")
    print("to the statement in the user's question, proving the answer is 'No'.")

if __name__ == '__main__':
    main()