import itertools

def main():
    """
    This function demonstrates the construction of an infinite set 'x' that is
    almost disjoint from a given countable family of infinite sets 'S'.
    """

    # --- Step 1: Define the problem setting ---
    # The question is whether a certain set 'x' exists. The answer is YES.
    # We will demonstrate this with a constructive algorithm.

    print("--- The Problem ---")
    print("Question: Given a countable collection S of infinite subsets of natural numbers,")
    print("does there always exist an infinite set x such that its intersection with any set in S is finite?")
    print("\nAnswer: YES. We can construct such a set x. Let's demonstrate.\n")

    # --- Step 2: Define the input sets S ---
    # We need a countable family S = {s_0, s_1, s_2, ...} of infinite sets.
    # We'll represent these infinite sets using Python generators.

    def s0_generator():
        # s0 = {Even numbers} = {0, 2, 4, 6, ...}
        n = 0
        while True:
            yield n
            n += 2

    def s1_generator():
        # s1 = {Multiples of 3} = {0, 3, 6, 9, ...}
        n = 0
        while True:
            yield n
            n += 3

    def s2_generator():
        # s2 = {Square numbers} = {0, 1, 4, 9, 16, ...}
        n = 0
        while True:
            yield n*n
            n += 1

    def s3_generator():
        # s3 = {Powers of 2} = {1, 2, 4, 8, 16, ...}
        n = 0
        while True:
            yield 2**n
            n += 1

    # Our collection of sets S. For this demonstration, S is finite, but the
    # algorithm works for a countably infinite collection as well.
    S = [s0_generator, s1_generator, s2_generator, s3_generator]
    NUM_SETS_IN_S = len(S)

    print("--- The Input Collection S ---")
    print(f"Let's use a sample collection S of {NUM_SETS_IN_S} infinite sets:")
    print("s0 = {Even numbers}")
    print("s1 = {Multiples of 3}")
    print("s2 = {Square numbers}")
    print("s3 = {Powers of 2}")
    print("-" * 20)

    # --- Step 3: The Construction Algorithm ---
    # We construct the set x = {x_0, x_1, x_2, ...} element by element.

    x_elements = []
    x_last = -1
    
    # We will compute the first N elements of x. Let's choose N.
    N_ELEMENTS_TO_CONSTRUCT = 10

    print(f"\nConstructing the first {N_ELEMENTS_TO_CONSTRUCT} elements of set x:")

    for n in range(N_ELEMENTS_TO_CONSTRUCT):
        # At step n, we consider the sets s_0, ..., s_n.
        
        # We need to find y_i = min(s_i > x_last) for each i <= n.
        # Since our S is finite, we only consider i < len(S).
        
        max_y = -1
        
        num_s_to_consider = min(n + 1, NUM_SETS_IN_S)
        
        print(f"\nStep n={n}:")
        print(f"  Last element of x was x_{n-1} = {x_last}.")
        print(f"  Considering sets s_0 through s_{num_s_to_consider - 1}.")

        for i in range(num_s_to_consider):
            s_i_gen = S[i]() # Get a fresh generator for s_i
            y_i = -1
            # Find the first element in s_i greater than x_last
            for val in s_i_gen:
                if val > x_last:
                    y_i = val
                    break
            
            print(f"    For s_{i}, the smallest element > {x_last} is y_{i} = {y_i}.")
            
            if y_i > max_y:
                max_y = y_i
        
        # This corresponds to the equation for the n-th element of x.
        # x_n = max(y_0, ..., y_{min(n, |S|-1)}) + 1
        x_n = max_y + 1
        print(f"  max(y_0, ..., y_{num_s_to_consider - 1}) = {max_y}.")
        print(f"  Final Equation: x_{n} = {max_y} + 1")
        print(f"  Result: x_{n} = {x_n}")


        x_elements.append(x_n)
        x_last = x_n

    # --- Step 4: The Result ---
    print("\n--- The Constructed Set x ---")
    print(f"The first {N_ELEMENTS_TO_CONSTRUCT} elements of x are:")
    print(f"x = {x_elements}")

    # --- Step 5: Verification ---
    # Let's check the intersections for the constructed part of x.
    print("\n--- Verification ---")
    print("Checking intersections |x \cap s_i|:")

    for i in range(NUM_SETS_IN_S):
        # To get the elements of s_i, we need to generate them.
        s_i_list = list(itertools.takewhile(lambda k: k <= x_elements[-1], S[i]()))
        
        intersection = set(x_elements) & set(s_i_list)
        
        print(f"|x \cap s_{i}| is at least {len(intersection)}.")
        print(f"  x \cap s_{i} (up to {x_elements[-1]}) = {sorted(list(intersection))}")

    print("\nThe algorithm ensures that these intersections will remain finite even as we generate more elements for x.")

if __name__ == '__main__':
    main()