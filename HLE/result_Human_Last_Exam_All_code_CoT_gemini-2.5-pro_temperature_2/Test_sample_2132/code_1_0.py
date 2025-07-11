import math

def calculate_boettcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Explain the chemistry and identify the product's formula.
    print("The Favorskii rearrangement of 2-chlorocyclohexanone yields the product methyl cyclopentanecarboxylate.")
    print("The chemical formula of this product is C7H12O2.")
    print("-" * 30)

    # Step 2: Define the atomic composition of the product molecule.
    atomic_composition = {'C': 7, 'H': 12, 'O': 2}
    n_C = atomic_composition['C']
    n_H = atomic_composition['H']
    n_O = atomic_composition['O']
    
    n_total = sum(atomic_composition.values())

    print("The Böttcher Molecular Complexity is calculated based on atomic composition.")
    print(f"Number of Carbon atoms (n_C) = {n_C}")
    print(f"Number of Hydrogen atoms (n_H) = {n_H}")
    print(f"Number of Oxygen atoms (n_O) = {n_O}")
    print(f"Total number of atoms (n_total) = {n_total}")
    print("-" * 30)

    # Step 3: Calculate the terms for the Böttcher formula.
    # C = n_total * log2(n_total) - sum(ni * log2(ni))
    term_total = n_total * math.log2(n_total)
    term_C = n_C * math.log2(n_C)
    term_H = n_H * math.log2(n_H)
    term_O = n_O * math.log2(n_O) if n_O > 0 else 0
    
    sum_terms = term_C + term_H + term_O
    complexity = term_total - sum_terms

    # Step 4: Display the calculation step-by-step.
    print("The formula is: C = n_total*log2(n_total) - (n_C*log2(n_C) + n_H*log2(n_H) + n_O*log2(n_O))")
    print("\nPlugging in the numbers for methyl cyclopentanecarboxylate:")
    
    equation = f"C = {n_total}*log2({n_total}) - ({n_C}*log2({n_C}) + {n_H}*log2({n_H}) + {n_O}*log2({n_O}))"
    print(equation)
    
    # Show intermediate values for clarity
    print(f"C = {term_total:.4f} - ({term_C:.4f} + {term_H:.4f} + {term_O:.4f})")
    print(f"C = {term_total:.4f} - {sum_terms:.4f}")

    # Step 5: Print the final result.
    print("-" * 30)
    print(f"The final Böttcher Molecular Complexity is: {complexity:.4f}")

if __name__ == "__main__":
    calculate_boettcher_complexity()