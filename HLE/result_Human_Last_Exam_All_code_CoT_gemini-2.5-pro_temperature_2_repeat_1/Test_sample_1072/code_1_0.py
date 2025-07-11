def calculate_he_age():
    """
    Calculates the (U-Th)/He ages for the three scenarios described.
    """
    # --- Constants ---
    T_surf = 25  # Surface temperature in °C
    geo_grad = 25  # Geothermal gradient in °C/km
    Tc_ZHe = 180   # Zircon (U-Th)/He closure temperature in °C
    Tc_AHe = 70    # Apatite (U-Th)/He closure temperature in °C

    # --- Sample 1: Zircon ---
    print("--- Sample 1: Zircon from Pluton ---")
    depth_s1_start = 15.0  # km
    time_s1_span = 100.0   # Ma
    # Depth at which temperature equals ZHe closure temperature
    depth_s1_closure = (Tc_ZHe - T_surf) / geo_grad
    # Assuming linear exhumation from depth_s1_start to 0 over time_s1_span
    # Age is proportional to depth: Age/time_span = Depth/depth_start
    age_s1 = (depth_s1_closure / depth_s1_start) * time_s1_span
    print(f"Initial depth: {depth_s1_start} km at {time_s1_span} Ma")
    print(f"Closure temperature (ZHe): {Tc_ZHe}°C")
    print(f"Depth of closure: ({Tc_ZHe} - {T_surf}) / {geo_grad} = {depth_s1_closure:.2f} km")
    print(f"Calculated Age = ({depth_s1_closure:.2f} / {depth_s1_start}) * {time_s1_span} = {age_s1:.2f} Ma\n")

    # --- Sample 2: Apatite ---
    print("--- Sample 2: Apatite from Sediment ---")
    T_s2_start = 250.0  # °C
    time_s2_span = 100.0  # Ma
    # Initial depth at 100 Ma
    depth_s2_start = (T_s2_start - T_surf) / geo_grad
    # Depth at which temperature equals AHe closure temperature
    depth_s2_closure = (Tc_AHe - T_surf) / geo_grad
    # Assuming linear exhumation from depth_s2_start to 0 over time_s2_span
    age_s2 = (depth_s2_closure / depth_s2_start) * time_s2_span
    print(f"Initial temperature: {T_s2_start}°C at {time_s2_span} Ma")
    print(f"Calculated initial depth: ({T_s2_start} - {T_surf}) / {geo_grad} = {depth_s2_start:.2f} km")
    print(f"Closure temperature (AHe): {Tc_AHe}°C")
    print(f"Depth of closure: ({Tc_AHe} - {T_surf}) / {geo_grad} = {depth_s2_closure:.2f} km")
    print(f"Calculated Age = ({depth_s2_closure:.2f} / {depth_s2_start:.2f}) * {time_s2_span} = {age_s2:.2f} Ma\n")
    
    # --- Sample 3: Apatite from Rhyolite ---
    print("--- Sample 3: Apatite from Rhyolite ---")
    age_s3 = 90.0 # Ma
    print(f"The sample cooled rapidly during eruption at {age_s3:.2f} Ma.")
    print(f"The age is the eruption age: {age_s3:.2f} Ma\n")
    
    # --- Final Ranking ---
    print("--- Age Ranking (Oldest to Youngest) ---")
    print(f"Sample 3 ({age_s3:.2f} Ma) > Sample 1 ({age_s1:.2f} Ma) > Sample 2 ({age_s2:.2f} Ma)")
    print("Therefore, statement [H] is correct.")

if __name__ == '__main__':
    calculate_he_age()