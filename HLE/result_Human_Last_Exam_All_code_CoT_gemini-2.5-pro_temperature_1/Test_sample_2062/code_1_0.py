import math

def solve_matrix_similarity_questions():
    """
    Solves the three-part question about the similarity of diagonal matrices.
    """
    # Part (a): Is it true that A and B are similar if and only if the
    # multiplicities of each eigenvalue are identical?
    # For diagonal matrices, this is true. A similarity class is uniquely
    # determined by the multiset of eigenvalues (diagonal entries).
    answer_a = "Yes"

    # Part (b): For n = 3 and distinct eigenvalues alpha, beta, gamma,
    # how many similarity classes exist?
    # Since the multiset of eigenvalues is fixed as {alpha, beta, gamma},
    # all such diagonal matrices belong to the same similarity class.
    answer_b = 1

    # Part (c): Does the number of similarity classes for diagonal matrices
    # in M_n(F) grow exponentially with n for fixed q?
    # The number of classes is the number of multisets of size n from a set
    # of size q, which is given by the binomial coefficient C(n + q - 1, n).
    # This is a polynomial in n of degree q-1, not an exponential function of n.
    answer_c = "No"

    # Print the final combined answer in the required format
    print("Final Answer:")
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

    # --- Demonstration for Part (c) ---
    # We compare the polynomial growth of the number of classes with exponential growth.
    print("\n" + "="*50)
    print("Demonstration for Part (c): Polynomial vs. Exponential Growth")
    print("="*50)
    
    # Let's fix q, the size of the field F_q.
    q = 4 # Example: a field with 4 elements.
    
    print(f"Number of similarity classes in M_n(F_q) for q={q}")
    print(f"Formula: C(n + {q-1}, n) = (n+{q-1})! / (n! * {q-1}!)")
    print("-" * 50)
    print(f"{'n':>3} | {'# Classes (Polynomial)':>25} | {'Exponential (e.g., 2^n)':>25}")
    print("-" * 50)

    for n in range(1, 16):
        # Polynomial growth: C(n + q - 1, n)
        # This is the number of similarity classes.
        num_classes = math.comb(n + q - 1, n)
        
        # Exponential growth for comparison
        exp_growth = 2**n
        
        print(f"{n:3d} | {num_classes:25,d} | {exp_growth:25,d}")

    print("-" * 50)
    print("As n increases, the exponential function grows much faster than the")
    print("polynomial function representing the number of similarity classes.")

# Execute the function to get the answers and demonstration.
solve_matrix_similarity_questions()