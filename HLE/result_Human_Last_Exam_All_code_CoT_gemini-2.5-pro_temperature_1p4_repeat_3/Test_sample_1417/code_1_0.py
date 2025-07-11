import pandas as pd

def find_optimal_material():
    """
    Compares materials based on their radiation length to find the optimal choice
    for a particle detector cooling system.
    
    The unique parameter for a particle detector is the need to minimize particle interaction,
    which is achieved by maximizing the radiation length (X0) of the material used.
    A longer radiation length means the material is more "transparent" to particles.
    """
    
    # Data for materials:
    # Z: Atomic Number
    # density: g/cm^3
    # X0_mass: Radiation Length in g/cm^2 (mass thickness)
    # The radiation length in cm (distance) is X0_mass / density.
    
    data = {
        'Titanium':         {'option': 'A', 'Z': 22, 'density': 4.51, 'X0_mass': 16.3},
        'Aluminium':        {'option': 'B', 'Z': 13, 'density': 2.70, 'X0_mass': 24.01},
        '316 Stainless Steel': {'option': 'C', 'Z': '~26', 'density': 8.00, 'X0_mass': 13.9}, # Using Iron as proxy
        'Copper':           {'option': 'D', 'Z': 29, 'density': 8.96, 'X0_mass': 12.86},
        'Nickel':           {'option': 'E', 'Z': 28, 'density': 8.90, 'X0_mass': 12.93},
    }

    # Use pandas for nice formatting
    df = pd.DataFrame(data).T
    
    # Calculate radiation length in cm
    df['X0_cm'] = df['X0_mass'] / df['density']

    print("--- Material Comparison for Particle Detector Cooling Tubes ---")
    print("The goal is to maximize the Radiation Length (X0) to minimize particle interaction.\n")
    
    # Sort by the radiation length in cm, descending
    df_sorted = df.sort_values(by='X0_cm', ascending=False)
    
    print(df_sorted[['option', 'Z', 'density', 'X0_cm']].round(2))

    # Identify the best choice
    best_material = df_sorted.index[0]
    best_option = df_sorted['option'][0]
    
    print(f"\nConclusion:")
    print(f"For pipes of similar dimensions, {best_material} is the optimum choice because it has the longest")
    print(f"radiation length ({df_sorted['X0_cm'][0]:.2f} cm), making it the most transparent to particles.")

find_optimal_material()
<<<B>>>