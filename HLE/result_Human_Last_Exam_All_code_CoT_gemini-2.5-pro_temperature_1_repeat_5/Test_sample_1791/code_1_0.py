import sys

def solve_electron_phonon_interaction():
    """
    This script derives the effective electron-electron interaction by performing
    a path integral over the phonon fields. It prints the derivation step-by-step.
    """
    # Use unicode for better rendering of mathematical symbols
    if sys.stdout.encoding.lower() not in ('utf-8', 'utf8'):
        print("Warning: Your terminal may not support Unicode characters. Formulas may not render correctly.")
        w, n, q, g, m, i, d = 'w', 'n', 'q', 'g', 'm', 'i', 'd'
        dagger = 't'
    else:
        w, n, q, g, m, i, d = '\u03c9', '\u03c1', 'q', 'g', 'm', 'i', '\u2202'
        dagger = '\u2020'

    print("Step 1: Define the Hamiltonian and Action")
    print("------------------------------------------")
    print(f"The phonon Hamiltonian is given by:")
    print(f"  H_ph = \u03A3_{{{q},j}} {w}_{q} a_{{{q},j}}^{dagger} a_{{{q},j}}")
    print("\nThe electron-phonon coupling is given as:")
    print(f"  H_e-ph = g \u03A3_{{k,{q},j}} ({i}{q}_j / (2{m}{w}_{q})\u00b9\u00b2) {n}_{q} (a_{{{q},j}} + a_{{-{q},j}}^{dagger})")
    print("\nNote: The provided H_e-ph is not Hermitian. We will proceed assuming the standard Hermitian form:")
    print(f"  H_e-ph = \u03A3_{{{q},j}} (C_{{{q}j}} {n}_{{-{q}}} a_{{{q}j}} + C_{{{q}j}}* {n}_{q} a_{{{q}j}}^{dagger})")
    print(f"where the coupling constant C_{{{q}j}} is:")
    print(f"  C_{{{q}j}} = g * ({i}{q}_j / \u221a(2{m}{w}_{q}))")
    print(f"and {n}_{q} is the electron density operator.")

    print("\nStep 2: Path Integral Formulation in Imaginary Time")
    print("--------------------------------------------------")
    print("We express the partition function Z as a path integral over electron and phonon fields.")
    print("The goal is to integrate out the phonon fields (a, a*) to find an effective action for the electrons.")
    print("In Matsubara frequency space ({w}_n), the action for the phonons and the interaction is:")
    print(f"  S_ph = \u03A3_{{{q},j,n}} a*_{{{q}j,n}} (-{i}{w}_n + {w}_{q}) a_{{{q}j,n}}")
    print(f"  S_int = \u03A3_{{{q},j,n}} [C_{{{q}j}} {n}_{{-{q},-n}} a_{{{q}j,n}} + C*_{{{q}j}} {n}_{{{q},n}} a*_{{{q}j,n}}]")

    print("\nStep 3: Gaussian Integration over Phonon Fields")
    print("------------------------------------------------")
    print("The total action involving phonons, S_ph + S_int, is quadratic in the phonon fields a and a*.")
    print("This allows us to perform the Gaussian integral exactly by 'completing the square'.")
    print("For each mode (q, j, n), we integrate \u222B da*da exp[-(A a*a + J*a + J a*)].")
    print("The result of this integration gives an effective action for the electron densities:")
    print(f"  S_eff = - (1/2) \u03A3_{{{q},j,n}} |C_{{{q}j}}|\u00b2 {n}_{{{q},n}}{n}_{{-{q},-n}} * (2{w}_{q} / ({w}_{q}\u00b2 + {w}_n\u00b2))")
    print("(The factor of 1/2 arises from carefully summing over q and -q modes).")

    print("\nStep 4: Identify the Effective Electron-Electron Interaction")
    print("-------------------------------------------------------------")
    print("The effective action S_eff can be written as:")
    print(f"  S_eff = (1/2) \u03A3_{{{q},n}} V_eff({q}, {i}{w}_n) {n}_{{{q},n}}{n}_{{-{q},-n}}")
    print("By comparing the two expressions for S_eff, we identify the effective interaction V_eff:")
    print(f"  V_eff({q}, {i}{w}_n) = - \u03A3_j |C_{{{q}j}}|\u00b2 * (2{w}_{q} / ({w}_{q}\u00b2 + {w}_n\u00b2))")

    print("\nStep 5: Substitute the Coupling Constant and Simplify")
    print("-------------------------------------------------------")
    print(f"First, we calculate the magnitude squared of the coupling constant C_{{{q}j}}:")
    print(f"  |C_{{{q}j}}|\u00b2 = (g * {i}{q}_j / \u221a(2{m}{w}_{q})) * (g * (-{i}{q}_j) / \u221a(2{m}{w}_{q}))")
    print(f"  |C_{{{q}j}}|\u00b2 = g\u00b2 * {q}_j\u00b2 / (2{m}{w}_{q})")
    print("\nNow, substitute this into the expression for V_eff:")
    print(f"  V_eff({q}, {i}{w}_n) = - \u03A3_j [g\u00b2 * {q}_j\u00b2 / (2{m}{w}_{q})] * [2{w}_{q} / ({w}_{q}\u00b2 + {w}_n\u00b2)]")
    print("The terms '2' and 'w_q' cancel out:")
    print(f"  V_eff({q}, {i}{w}_n) = - \u03A3_j [g\u00b2 * {q}_j\u00b2 / m] * [1 / ({w}_{q}\u00b2 + {w}_n\u00b2)]")
    print("\nAssuming the phonon frequency {w}_{q} is independent of the branch j and summing over the spatial components j (e.g., j=x,y,z),")
    print(f"we get \u03A3_j {q}_j\u00b2 = |{q}|\u00b2.")

    print("\nFinal Result")
    print("------------")
    print("The effective electron-electron interaction for a given momentum q is:")
    final_result = f"V_eff(q, i{w}_n) = - (g\u00b2 |q|\u00b2) / (m * ({w}_{q}\u00b2 + {w}_n\u00b2))"
    print(f"  {final_result}")
    print("\nThis is an attractive interaction (V_eff < 0) that depends on the momentum transfer q")
    print("and the Matsubara frequency {w}_n, which corresponds to the energy transfer.")

if __name__ == '__main__':
    solve_electron_phonon_interaction()
