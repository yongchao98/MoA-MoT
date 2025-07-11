def display_bare_greens_function():
    """
    This function prints the symbolic formula for the bare Green's function G_0,
    showing its dependence on the single-particle energy epsilon_k.
    """
    
    # Define the symbols in the equation for clarity.
    # These are not computational variables but labels for the printed output.
    numerator = "1"
    real_frequency = "ω"
    single_particle_energy = "ϵ_k"
    chemical_potential = "μ"
    infinitesimal = "iδ" # 'i' is the imaginary unit, 'δ' is a positive infinitesimal
    
    # Construct the components of the equation string
    G_0 = f"G_0(k, {real_frequency})"
    denominator = f"{real_frequency} - ({single_particle_energy} - {chemical_potential}) + {infinitesimal}"
    
    # The standard form for the time-ordered Green's function is often written as:
    # G_0(k, ω) = 1 / (ω - ϵ_k + iδ sgn(ϵ_k - μ))
    # A more common and general form for the retarded Green's function is shown below.
    # For simplicity, we will display the general relationship which captures the core dependence.
    denominator_standard = f"({real_frequency} - {single_particle_energy} + {infinitesimal})"
    
    # To be precise, let's represent the dependence for the retarded Green's function,
    # obtained via analytic continuation iω_n -> ω + iδ, which results in:
    # G_0(k, ω) = 1 / ( (ϵ_k - μ) - ω - iδ ) which is equivalent to 1 / ( ω - (ϵ_k - μ) + iδ' )
    # Let's print the most common form found in textbooks.
    final_denominator = f"{real_frequency} - {single_particle_energy} + {infinitesimal.replace('δ', 'δ * sgn(ϵ_k - μ)')}"
    # However, for showing the fundamental dependence on ϵ_k, the simplest form is clearest.
    
    print("The bare Green's function G_0 has the following dependence on the single-particle energy ϵ_k:")
    print("\n")
    print(f"         {numerator}")
    print(f"{G_0} = {'-' * (len(denominator) + 2)}")
    print(f"         {denominator}")
    print("\n")
    print("Where:")
    print(f"  G_0: The bare Green's function for momentum k and frequency {real_frequency}")
    print(f"  {real_frequency}: The frequency (energy) of the particle")
    print(f"  {single_particle_energy}: The single-particle energy eigenvalue")
    print(f"  {chemical_potential}: The chemical potential")
    print(f"  {infinitesimal}: An infinitesimal term that defines the pole structure")

display_bare_greens_function()