import sys

def solve_fermi_hubbard_properties():
    """
    This function outlines the properties of a 1D Fermi-Hubbard system
    with on-site two-body losses in the infinite time limit, based on physical reasoning.
    """

    # The reasoning leads to a final state that is a superposition of "jammed"
    # configurations. We analyze the properties of such a state.

    # 1) Zero tunneling: TRUE. Jammed states are eigenstates of the tunneling
    #    operator with eigenvalue zero, as hopping is Pauli-blocked.
    prop1 = True

    # 2) Zero particles: FALSE. The system can get stuck in a non-vacuum jammed
    #    state before all particles are lost.
    prop2 = False

    # 3) Zero losses: TRUE. The final state must be a steady state, which by
    #    definition has no ongoing particle loss. This implies zero double occupancy.
    prop3 = True

    # 4) Spin entanglement: TRUE. The final state can be a superposition of
    #    degenerate, jammed ground states (e.g., |udud> and |dudu>),
    #    which results in an entangled state.
    prop4 = True

    # 5) Anti-ferromagnetic-like spin correlations: TRUE. The dissipative loss of
    #    high-energy doublons acts like cooling, driving the system to its
    #    low-energy state, which exhibits AF correlations.
    prop5 = True

    # 6) Ferromagnetic-like spin correlations: FALSE. These are energetically
    #    unfavorable compared to anti-ferromagnetic correlations.
    prop6 = False

    properties = {
        1: ("Zero tunneling", prop1),
        2: ("Zero particles", prop2),
        3: ("Zero losses", prop3),
        4: ("Spin entanglement", prop4),
        5: ("Anti-ferromagnetic-like spin correlations", prop5),
        6: ("Ferromagnetic-like spin correlations", prop6),
    }

    final_property_indices = [k for k, v in properties.items() if v[1]]

    print("The final state has the following properties:")
    for index in final_property_indices:
        print(f"{index}) {properties[index][0]}")

    print("\nThe final combination of properties is represented by the equation:")
    # Formatting as an "equation" as requested by the prompt
    equation_str = " + ".join(map(str, final_property_indices))
    print(equation_str)


if __name__ == "__main__":
    solve_fermi_hubbard_properties()
