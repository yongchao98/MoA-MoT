import textwrap

def analyze_boronic_acid_synthesis():
    """
    Analyzes a problematic boronic acid synthesis and suggests a solution.
    """

    # --- Problem Statement ---
    reaction_description = """
    A reaction was performed to synthesize (2-bromo-4-chlorophenyl)boronic acid.
    The issue is that the final product shows two different boron (11B NMR) signals,
    whereas a single pure product should only show one.
    """
    print("--- The Problem ---")
    print(textwrap.dedent(reaction_description).strip())

    # --- Analysis of the Reaction Stoichiometry ---
    analysis = """
    Let's analyze the reaction equation and stoichiometry used:

    Problematic Reaction:
    1 eq (2-bromo-4-chloro-1-iodobenzene) + 1.05 eq (n-BuLi) + 5 eq (trimethyl borate)
    -->  Mixture of products (giving two 11B NMR signals)

    The key detail here is the use of 5 equivalents of trimethyl borate. This is a very large excess.
    The typical amount used in these reactions is about 1.2 to 2.0 equivalents.

    This large excess is problematic because upon aqueous workup, all unused trimethyl borate is
    converted into boric acid (B(OH)3). The final material is therefore a mixture of the
    desired (2-bromo-4-chlorophenyl)boronic acid and a large amount of boric acid. These
    two different boron-containing compounds will result in two separate signals in the 11B NMR.
    """
    print("\n--- Stoichiometry Analysis ---")
    print(textwrap.dedent(analysis).strip())

    # --- Proposed Solution ---
    solution = """
    The most direct way to solve this problem is to reduce the amount of the electrophile,
    trimethyl borate, to a more standard and reasonable excess. This minimizes the formation
    of boric acid byproduct, making purification easier and leading to a single, clean product.

    Corrected Reaction:
    1 eq (2-bromo-4-chloro-1-iodobenzene) + 1.05 eq (n-BuLi) + 1.5 eq (trimethyl borate)
    -->  (2-bromo-4-chlorophenyl)boronic acid (single product, single 11B NMR signal)

    Therefore, the best course of action is to use less trimethyl borate.
    """
    print("\n--- Proposed Solution ---")
    print(textwrap.dedent(solution).strip())


analyze_boronic_acid_synthesis()
<<<D>>>