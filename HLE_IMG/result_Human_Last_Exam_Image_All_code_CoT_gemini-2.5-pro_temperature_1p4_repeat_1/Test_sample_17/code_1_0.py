import math

def calculate_theoretical_stress():
    """
    This function explains and calculates the theoretical stress at the tip of a
    perfectly sharp notch based on the principles of linear elasticity.
    """
    
    # Define variables as strings for the explanation
    nominal_stress = "σ_y"
    max_stress = "σ_max"
    stress_concentration_factor = "K_t"
    tip_radius = "ρ"

    print("The problem is to find the theoretical stress at the tip of a sharp wedge (Point A).")
    print("This is a classic problem of stress concentration in linear elasticity.\n")

    print(f"1. The relationship between maximum stress ({max_stress}) at a notch and the nominal applied stress ({nominal_stress}) is given by:")
    print(f"   {max_stress} = {stress_concentration_factor} * {nominal_stress}\n")

    print(f"2. The stress concentration factor ({stress_concentration_factor}) for a notch is highly dependent on the sharpness of its tip, specifically the tip radius ({tip_radius}).\n")

    print(f"3. For a theoretically perfect, sharp wedge as shown in the diagram, the tip radius ({tip_radius}) is considered to be zero.")
    print(f"   {tip_radius} = 0\n")

    print(f"4. As the tip radius ({tip_radius}) approaches zero, the stress concentration factor ({stress_concentration_factor}) approaches infinity.")
    print(f"   {stress_concentration_factor} -> ∞\n")
    
    print("5. Therefore, we can write the final equation for the theoretical stress at the tip:")
    
    # We use the concept of infinity to represent the theoretical value
    final_kt = "∞"
    
    # The 'equation' part requested by the user
    print(f"   {max_stress} = {final_kt} * {nominal_stress}\n")
    
    print("Conclusion: The theoretical stress at the tip of a perfectly sharp wedge is infinite.")

if __name__ == '__main__':
    calculate_theoretical_stress()