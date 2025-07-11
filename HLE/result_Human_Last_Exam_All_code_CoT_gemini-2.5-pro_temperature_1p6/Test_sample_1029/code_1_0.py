import math

def display_poynting_vector():
    """
    This function prints the derived symbolic formulas for the Poynting vector
    both inside and outside the described cylindrical rod.
    """
    
    # Symbolic placeholders for the variables
    rho = "ρ"
    E = "E"
    v = "v"
    R = "R"
    r = "r"
    epsilon_0 = "ε₀"
    
    print("The Poynting vector S has two components: a radial component (S_r) and an axial component (S_z).")
    print("S = (S_r) r_hat + (S_z) z_hat\n")

    # --- Inside the rod (r < R) ---
    print("Inside the rod (where r < R):")
    
    # Radial component S_r
    # The numbers in the equation are 1 (implicit) and 2.
    S_r_in_num_1 = E
    S_r_in_num_2 = rho
    S_r_in_num_3 = v
    S_r_in_num_4 = r
    S_r_in_den_1 = 2
    print(f"S_r = -({S_r_in_num_1} * {S_r_in_num_2} * {S_r_in_num_3} * {S_r_in_num_4}) / {S_r_in_den_1}")
    
    # Axial component S_z
    # The numbers in the equation are 2, 2, 4.
    S_z_in_num_1 = rho
    S_z_in_num_power_1 = 2
    S_z_in_num_2 = v
    S_z_in_num_3 = r
    S_z_in_num_power_2 = 2
    S_z_in_den_1 = 4
    S_z_in_den_2 = epsilon_0
    print(f"S_z = ({S_z_in_num_1}**{S_z_in_num_power_1} * {S_z_in_num_2} * {S_z_in_num_3}**{S_z_in_num_power_2}) / ({S_z_in_den_1} * {S_z_in_den_2})")
    
    print("-" * 50)
    
    # --- Outside the rod (r > R) ---
    print("Outside the rod (where r > R):")

    # Radial component S_r
    # The numbers in the equation are 2, 2.
    S_r_out_num_1 = E
    S_r_out_num_2 = rho
    S_r_out_num_3 = v
    S_r_out_num_4 = R
    S_r_out_num_power_1 = 2
    S_r_out_den_1 = 2
    S_r_out_den_2 = r
    print(f"S_r = -({S_r_out_num_1} * {S_r_out_num_2} * {S_r_out_num_3} * {S_r_out_num_4}**{S_r_out_num_power_1}) / ({S_r_out_den_1} * {S_r_out_den_2})")

    # Axial component S_z
    # The numbers in the equation are 2, 4, 4, 2.
    S_z_out_num_1 = rho
    S_z_out_num_power_1 = 2
    S_z_out_num_2 = v
    S_z_out_num_3 = R
    S_z_out_num_power_2 = 4
    S_z_out_den_1 = 4
    S_z_out_den_2 = epsilon_0
    S_z_out_den_3 = r
    S_z_out_den_power_1 = 2
    print(f"S_z = ({S_z_out_num_1}**{S_z_out_num_power_1} * {S_z_out_num_2} * {S_z_out_num_3}**{S_z_out_num_power_2}) / ({S_z_out_den_1} * {S_z_out_den_2} * {S_z_out_den_3}**{S_z_out_den_power_1})")

display_poynting_vector()