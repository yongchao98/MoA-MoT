def solve_ccsd_question():
    """
    Explains for which excited Slater determinants the CCSD projection equation is zero.
    
    The problem concerns the matrix element <Φ_I|H_bar|Φ> in CCSD theory.
    
    H_bar = exp(-T) * H * exp(T)
    T = T1 + T2 (for CCSD)
    H contains at most 2-body interactions.
    
    The amplitudes are solved by enforcing <Φ_I|H_bar|Φ> = 0 for single (S)
    and double (D) excitations. The question is for which OTHER excitations this
    is also zero.
    """

    print("Step 1: Understanding the operator H_bar")
    print("The similarity-transformed Hamiltonian, H_bar, can be expressed as a sum of connected diagrams involving one H vertex and any number of T1 and T2 vertices.")
    print("H_bar = (H * exp(T))_connected")
    print("-" * 30)

    print("Step 2: Analyzing the 'body rank' of H_bar")
    print("The Hamiltonian H is a two-body operator. The T1 and T2 operators describe one- and two-body correlations.")
    print("When combined in connected diagrams, they create 'effective' many-body interactions.")
    print("A two-body H can connect to T2 operators to create higher-order interactions.")
    print("\n- An H operator alone is a 2-body interaction, creating up to double excitations.")
    print("- A connected H-T2 term creates an effective 3-body interaction, creating up to triple excitations.")
    print("- A connected H-T2-T2 term creates an effective 4-body interaction, creating up to quadruple excitations.")
    print("\nBecause the fundamental interaction H is two-body, it cannot form a connected diagram that represents a five-body or higher interaction.")
    print("Thus, the maximum body rank of H_bar in CCSD is 4.")
    print("-" * 30)

    print("Step 3: Calculating the Maximum Excitation Level")
    print("An n-body operator can connect a reference determinant to determinants that are at most n-tuply excited.")
    max_body_rank_H = 2
    # Each connection to an independent T2 operator can increase the effective body rank by one.
    # An H vertex can link at most two independent T2 clusters.
    additional_rank_from_T2_a = 1
    additional_rank_from_T2_b = 1
    
    max_excitation_level = max_body_rank_H + additional_rank_from_T2_a + additional_rank_from_T2_b
    
    print(f"The maximum excitation level M for which <Φ_M|H_bar|Φ> can be non-zero is given by:")
    print(f"M_max = (Body Rank of H) + (Max additional rank from correlations)")
    print(f"M_max = {max_body_rank_H} + {additional_rank_from_T2_a} + {additional_rank_from_T2_b} = {max_excitation_level}")
    print("-" * 30)
    
    print("Step 4: Conclusion")
    print(f"Since H_bar is at most a {max_excitation_level}-body operator, it cannot connect the reference determinant |Φ> to any determinant |Φ_I> that is excited by more than {max_excitation_level} orbitals.")
    print(f"\nTherefore, the matrix element <Φ_I|H_bar|Φ> is guaranteed to be zero for:")
    print("All Slater determinants corresponding to Quintuple (5), Sextuple (6), and all higher-level excitations.")

solve_ccsd_question()
