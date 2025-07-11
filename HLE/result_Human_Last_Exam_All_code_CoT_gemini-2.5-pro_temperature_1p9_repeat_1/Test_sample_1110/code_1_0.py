import textwrap

def solve_chemistry_problem():
    """
    This script analyzes a chemical reaction problem and provides a step-by-step solution.
    """
    print("Analyzing the chemical problem step-by-step:")
    print("-" * 50)

    # Step 1: Define the problem
    problem_description = """
    The reaction of 2-bromo-4-chloro-1-iodobenzene with 1.05 equivalents of n-BuLi and
    5 equivalents of trimethyl borate was performed. The NMR analysis showed two
    distinct boron signals, indicating the formation of two different boronic acid products
    instead of just the single intended product.
    """
    print("1. The Problem:")
    print(textwrap.indent(textwrap.dedent(problem_description), '    '))

    # Step 2: Analyze the reaction pathway
    analysis = """
    The intended reaction is a selective lithium-halogen exchange at the most reactive
    site, which is the carbon-iodine (C-I) bond. The established reactivity order
    for this exchange is I > Br > Cl.

    Intended Reaction:
    2-bromo-4-chloro-1-iodobenzene + 1.0 eq n-BuLi  -> (2-bromo-4-chlorophenyl)lithium

    However, the presence of a second product suggests a side reaction is occurring.
    The most plausible side reaction involves the reaction of n-BuLi at the second-most
    reactive site, the carbon-bromine (C-Br) bond.
    """
    print("\n2. Reaction Analysis:")
    print(textwrap.indent(textwrap.dedent(analysis), '    '))

    # Step 3: Identify the cause
    cause = """
    This lack of selectivity is almost certainly caused by the stoichiometry of n-BuLi.
    The reaction uses 1.05 equivalents, which is a slight excess. This excess n-BuLi can
    react with the C-Br bond after the more reactive C-I bond has been fully consumed.
    Inaccuracies in determining n-BuLi concentration, a common issue, can worsen this
    problem. The excess trimethyl borate (5 eq) is a trapping agent and is not the
    cause of the formation of the second organolithium species.
    """
    print("\n3. Root Cause Identification:")
    print(textwrap.indent(textwrap.dedent(cause), '    '))

    # Step 4: Propose the solution
    solution = """
    To prevent the side reaction, the amount of n-BuLi must be carefully controlled.
    Using a precisely measured, stoichiometric amount (exactly 1.00 eq, confirmed by
    recent titration) will ensure that the n-BuLi is consumed by the more reactive C-I
    bond, leaving none to react with the less reactive C-Br bond. This will prevent the
    formation of the undesired side product.
    
    The proposed ideal equation would be:
    - Step 1: 2-bromo-4-chloro-1-iodobenzene + 1.00 eq n-BuLi -> ...
    - Step 2: ... + >1.0 eq trimethyl borate -> Final Product
    """
    print("\n4. Conclusion and Solution:")
    print(textwrap.indent(textwrap.dedent(solution), '    '))

    final_answer = 'C'
    print(f"\nBased on the analysis, the best solution is to use a more precise amount of n-BuLi.")

if __name__ == '__main__':
    solve_chemistry_problem()
<<<C>>>