import itertools

def calculate_b3_variance():
    """
    Calculates the variance of the Coxeter length random variable on the 
    hyperoctahedral group of rank 3 (B_3).
    """
    n = 3

    # Generate permutations of {1, 2, ..., n}
    base_permutations = list(itertools.permutations(range(1, n + 1)))

    # Generate all signed permutations for B_n
    # An element w is represented as a tuple (w(1), w(2), ...)
    signed_permutations = []
    for p in base_permutations:
        for i in range(2**n):
            signs = []
            temp_i = i
            # Generate the tuple of n signs, e.g. (1, 1, -1)
            for _ in range(n):
                if temp_i % 2 == 0:
                    signs.append(1)
                else:
                    signs.append(-1)
                temp_i //= 2
            
            element = tuple(p[j] * signs[j] for j in range(n))
            signed_permutations.append(element)

    def coxeter_length(w):
        """
        Calculates the Coxeter length of a signed permutation w in B_n.
        l(w) = inv'(w) + nsp(w) + neg(w)
        """
        size = len(w)
        inv_prime = 0
        nsp = 0
        for i in range(size):
            for j in range(i + 1, size):
                if w[i] > w[j]:
                    inv_prime += 1
                if w[i] + w[j] < 0:
                    nsp += 1
        
        neg = sum(1 for x in w if x < 0)
        
        return inv_prime + nsp + neg

    # Calculate lengths for all elements
    all_lengths = [coxeter_length(w) for w in signed_permutations]

    # Calculate the components for the variance formula
    num_elements = len(all_lengths)
    sum_of_lengths = sum(all_lengths)
    mean_length = sum_of_lengths / num_elements

    sum_of_squares = sum(l**2 for l in all_lengths)
    mean_of_squares = sum_of_squares / num_elements
    
    variance = mean_of_squares - mean_length**2
    
    # Print the values required for the final equation
    print("To find the variance, we compute E[L^2] - (E[L])^2.")
    print(f"The number of elements |B_3| is {num_elements}.")
    print(f"The sum of lengths ΣL is {sum_of_lengths}.")
    print(f"The mean length E[L] is {sum_of_lengths} / {num_elements} = {mean_length}.")
    print(f"The sum of squared lengths Σ(L^2) is {sum_of_squares}.")
    print(f"The mean of squared lengths E[L^2] is {sum_of_squares} / {num_elements} = {mean_of_squares}.")
    print("\nThe final equation for the variance is:")
    print(f"{mean_of_squares} - ({mean_length})^2 = {variance}")


calculate_b3_variance()