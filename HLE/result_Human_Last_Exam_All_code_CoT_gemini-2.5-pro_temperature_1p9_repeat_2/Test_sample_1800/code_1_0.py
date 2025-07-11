import textwrap

def get_optimal_ni_ce_ratio_info():
    """
    This function provides information on the ideal Ni/Ce ratio for catalysis
    based on scientific literature.
    """
    
    # --- Ideal Ratio Based on Literature ---
    # The ideal ratio is often expressed as an atomic ratio. Studies show an
    # optimal range, with a Ni:Ce ratio of 1:4 being highly effective for
    # reactions like the Water Gas Shift (WGS).
    ni_atoms = 1
    ce_atoms = 4
    
    # --- Explanation ---
    explanation = f"""
The ideal atomic ratio of Nickel (Ni) to Cerium (Ce) for Ni-Ceria catalysts is crucial for maximizing performance in reactions like the Water Gas Shift (WGS) and water splitting. While the exact optimum can vary, a widely cited effective ratio is {ni_atoms}:{ce_atoms}.

Scientific Rationale:
1.  Maximizing Active Interface: The catalytic activity heavily relies on the interface between Ni nanoparticles and the Ceria (CeO₂) support. A ratio of {ni_atoms}:{ce_atoms} promotes high dispersion of Ni, creating a large interfacial area.
2.  Preventing Sintering: Higher concentrations of Nickel can lead to the agglomeration (sintering) of Ni particles at high reaction temperatures, which deactivates the catalyst by reducing the number of active sites. A lower Ni content helps maintain small, stable nanoparticles.
3.  Synergistic Effect: This ratio facilitates the synergistic effect where the Ceria support provides oxygen vacancies that participate in the catalytic cycle, and the Ni provides sites for reactant adsorption and dissociation.

Therefore, an equation representing this ideal ratio is:
"""
    
    # --- Print the results ---
    print(textwrap.dedent(explanation).strip())
    
    ratio = ni_atoms / ce_atoms
    print("\nEquation for the ideal atomic ratio:")
    print(f"Ni atoms / Ce atoms = {ni_atoms} / {ce_atoms} = {ratio}")

if __name__ == '__main__':
    get_optimal_ni_ce_ratio_info()