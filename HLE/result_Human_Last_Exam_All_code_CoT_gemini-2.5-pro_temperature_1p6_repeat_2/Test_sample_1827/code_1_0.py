import math

def calculate_goldstone_bosons():
    """
    Calculates the number of Goldstone bosons for a QCD phase transition.
    
    The problem considers a system of Nf quarks where one quark (strange) has a
    chemical potential, leading to a phase transition. We calculate the number of
    broken generators to find the number of Goldstone bosons.
    """
    
    # For a K-meson, the relevant quarks are up, down, and strange, so Nf = 3.
    Nf = 3

    print("--- Goldstone Boson Calculation in QCD ---")
    print(f"The system is described with an initial number of quark flavors, Nf = {Nf}.")
    print("-" * 40)

    # --- Step 1: Symmetry in the Gas Phase ---
    print("\n[Phase 1: Gas Phase (Before Condensation)]")
    
    # In the gas phase, the Lagrangian has a chemical potential for one quark (strange quark).
    # This, along with its different mass, breaks the SU(Nf) flavor symmetry.
    # The remaining ("iso-vector") symmetry is the subgroup that does not mix the strange
    # quark with the other (Nf-1) quarks. This group is G = SU(Nf-1) x U(1).
    # For Nf=3, this is G = SU(2) x U(1).
    
    print("The symmetry group G of the gas phase is SU(Nf-1) x U(1).")
    print(f"For Nf={Nf}, this is G = SU({Nf-1}) x U(1).")

    # Calculate the number of generators for G.
    # dim(SU(N)) = N^2 - 1
    # dim(U(1)) = 1
    dim_su_nf_minus_1 = (Nf - 1)**2 - 1
    dim_u1 = 1
    dim_G = dim_su_nf_minus_1 + dim_u1
    
    print("\nCalculating the number of generators for group G:")
    print(f"Number of generators for SU({Nf-1}) = ({Nf-1})^2 - 1 = {dim_su_nf_minus_1}")
    print(f"Number of generators for U(1) = {dim_u1}")
    print(f"Total generators for G = {dim_su_nf_minus_1} + {dim_u1} = {dim_G}")
    print("-" * 40)

    # --- Step 2: Symmetry in the Condensed Phase ---
    print("\n[Phase 2: Condensed Phase (After Condensation)]")
    
    # The problem states that the strange quark condenses, and the system effectively
    # behaves like one with (Nf-1) quarks.
    # The symmetry of this new effective system is the flavor symmetry of the remaining quarks, H = SU(Nf-1).
    # For Nf=3, this is H = SU(2).
    
    print("The symmetry group H of the condensed phase is SU(Nf-1).")
    print(f"For Nf={Nf}, this is H = SU({Nf-1}).")

    # Calculate the number of generators for H.
    dim_H = (Nf - 1)**2 - 1
    
    print("\nCalculating the number of generators for group H:")
    print(f"Number of generators for SU({Nf-1}) = ({Nf-1})^2 - 1 = {dim_H}")
    print("-" * 40)

    # --- Step 3: Calculate Goldstone Bosons ---
    print("\n[Result: Number of Goldstone Bosons]")
    
    # The number of Goldstone bosons is the number of broken generators, which is dim(G) - dim(H).
    num_goldstones = dim_G - dim_H
    
    print("According to Goldstone's theorem, the number of Goldstone bosons equals")
    print("the number of spontaneously broken symmetry generators.")
    print("\nFinal Equation:")
    print(f"Number of Goldstone Bosons = (Generators of G) - (Generators of H)")
    print(f"                           = {dim_G} - {dim_H} = {num_goldstones}")


if __name__ == '__main__':
    calculate_goldstone_bosons()