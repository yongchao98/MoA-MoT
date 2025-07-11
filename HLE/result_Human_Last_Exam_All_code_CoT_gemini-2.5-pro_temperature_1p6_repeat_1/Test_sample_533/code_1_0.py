def solve_h2_fock_space_decomposition():
    """
    This script explains the decomposition of the H2 molecule's Fock space
    in a minimal basis into symmetry-adapted Hilbert spaces and calculates
    the total number of such spaces.
    """
    print("Step-by-Step Decomposition of the H2 Fock Space")
    print("=================================================")

    print("\n1. System and Basis:")
    print("   - System: H2 molecule in a minimal basis (1s orbital per atom).")
    print("   - This creates 2 spatial Molecular Orbitals (MOs): σ_g (gerade) and σ_u (ungerade).")
    print("   - With spin, we have 4 spin-orbitals.")
    print("   - The full Fock space contains 2^4 = 16 possible electronic states.")

    print("\n2. Symmetries and Quantum Numbers:")
    print("   - The Hamiltonian commutes with operators for:")
    print("     - N (Total electron number)")
    print("     - S² (Total spin squared)")
    print("     - Spatial symmetry (D∞h point group)")
    print("   - A symmetry-adapted Hilbert space is a subspace where all states share")
    print("     a unique set of quantum numbers {N, S, Γ}, where Γ is the term symbol")
    print("     like ¹Σ_g⁺, which combines spin multiplicity (2S+1) and spatial symmetry.")

    print("\n3. Decomposition into N-electron Sectors:")
    
    total_spaces = 0
    sector_counts = []

    # N=0 Sector
    n0_spaces = 1
    total_spaces += n0_spaces
    sector_counts.append(n0_spaces)
    print("\n   --- Sector N = 0 (1 state) ---")
    print("   - Configuration: Vacuum (no electrons).")
    print("   - Term Symbol: ¹Σ_g⁺")
    print(f"   - Number of unique Hilbert spaces in this sector: {n0_spaces}")

    # N=1 Sector
    n1_spaces = 2
    total_spaces += n1_spaces
    sector_counts.append(n1_spaces)
    print("\n   --- Sector N = 1 (4 states) ---")
    print("   - Configurations: (σ_g)¹ or (σ_u)¹.")
    print("   - These result in two distinct term symbols: ²Σ_g⁺ and ²Σ_u⁺.")
    print(f"   - Number of unique Hilbert spaces in this sector: {n1_spaces}")
    
    # N=2 Sector
    n2_spaces = 3
    total_spaces += n2_spaces
    sector_counts.append(n2_spaces)
    print("\n   --- Sector N = 2 (6 states, for neutral H2) ---")
    print("   - Configurations: (σ_g)², (σ_u)², and (σ_g)¹(σ_u)¹.")
    print("   - (σ_g)² -> ¹Σ_g⁺")
    print("   - (σ_u)² -> ¹Σ_g⁺")
    print("   - (σ_g)¹(σ_u)¹ -> ¹Σ_u⁺ and ³Σ_u⁺")
    print("   - The unique symmetry spaces are defined by {S, Γ}:")
    print("     - ¹Σ_g⁺ (This space contains two configurations, making it 2D)")
    print("     - ¹Σ_u⁺")
    print("     - ³Σ_u⁺")
    print(f"   - Number of unique Hilbert spaces in this sector: {n2_spaces}")
    
    # N=3 Sector
    n3_spaces = 2
    total_spaces += n3_spaces
    sector_counts.append(n3_spaces)
    print("\n   --- Sector N = 3 (4 states) ---")
    print("   - This sector is related to N=1 by particle-hole symmetry.")
    print("   - Configurations correspond to one 'hole' in a full shell.")
    print("   - Term Symbols: ²Σ_g⁺ and ²Σ_u⁺.")
    print(f"   - Number of unique Hilbert spaces in this sector: {n3_spaces}")

    # N=4 Sector
    n4_spaces = 1
    total_spaces += n4_spaces
    sector_counts.append(n4_spaces)
    print("\n   --- Sector N = 4 (1 state) ---")
    print("   - Configuration: (σ_g)²(σ_u)² (fully occupied).")
    print("   - Term Symbol: ¹Σ_g⁺.")
    print(f"   - Number of unique Hilbert spaces in this sector: {n4_spaces}")

    print("\n4. Final Calculation:")
    print("   The maximum number of symmetry-adapted Hilbert spaces is the sum")
    print("   of the counts from each N-electron sector.")
    
    equation_str = " + ".join(map(str, sector_counts))
    print(f"   Equation: {sector_counts[0]} (for N=0) + {sector_counts[1]} (for N=1) + {sector_counts[2]} (for N=2) + {sector_counts[3]} (for N=3) + {sector_counts[4]} (for N=4) = {total_spaces}")
    
    print(f"\n   The maximum number of symmetry-adapted Hilbert spaces is {total_spaces}.")

solve_h2_fock_space_decomposition()

print("\n<<<9>>>")