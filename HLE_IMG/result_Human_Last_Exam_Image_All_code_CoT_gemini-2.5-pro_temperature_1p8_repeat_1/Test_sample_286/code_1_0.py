import math

def check_plot(plot_name, data):
    """
    Checks if a plot represents a physically valid quantum evolution based on estimated data.
    """
    print(f"--- Analyzing Plot {plot_name} ---")
    
    # Handle obvious violations first
    if 'invalid_values' in data:
        for point in data['invalid_values']:
            if point['name'] == 'sz' and abs(point['value']) > 1:
                print(f"Reason: Unphysical. Expectation value of σz cannot have magnitude > 1. Found <σz> ≈ {point['value']} at t≈{point['t']}.")
                return False
            if point['name'] == 'S' and point['value'] < 0:
                print(f"Reason: Unphysical. Entropy S cannot be negative. Found S ≈ {point['value']} at t≈{point['t']}.")
                return False

    # Extract initial conditions
    sz0 = data['t0']['sz']
    sp_mag0 = data['t0']['sp_mag']
    s0 = data['t0']['S']

    # The core physical condition is |r|² = <σz>² + 4 * |<σ+>|² <= 1
    # At t=0, entropy S=0, so the state must be pure, |r|²=1.
    print(f"At t=0, the plot shows S=0, which means the state must be pure, so |r|² should be 1.")
    print("Checking this condition: |r(0)|² = <σz>(0)² + 4 * |<σ+>|(0)²")
    r_sq_0 = sz0**2 + 4 * sp_mag0**2
    
    # We output each number in the final equation as requested
    print(f"Calculation: ({sz0})² + 4 * ({sp_mag0})² = {sz0**2:.2f} + 4 * {sp_mag0**2:.2f} = {sz0**2:.2f} + {4*sp_mag0**2:.2f} = {r_sq_0:.2f}")

    if not math.isclose(r_sq_0, 1.0, abs_tol=0.05): # Use a small tolerance for visual estimation
        print(f"Reason: Unphysical. For S=0, |r|² must be 1, but it is {r_sq_0:.2f}.")
        return False
        
    # Check that for t>0, the state remains inside the Bloch sphere.
    print("Condition at t=0 is met. Now checking for t > 0...")
    for t, point_data in data.get('t_other', []):
        sz_t = point_data.get('sz')
        sp_mag_t = point_data.get('sp_mag')
        r_sq_t = sz_t**2 + 4 * sp_mag_t**2
        print(f"At t={t}, |r|² = ({sz_t})² + 4 * ({sp_mag_t})² = {sz_t**2:.2f} + 4 * {sp_mag_t**2:.3f} = {sz_t**2:.2f} + {4*sp_mag_t**2:.3f} = {r_sq_t:.3f}")
        if r_sq_t > 1.0 + 0.05: # Add tolerance
            print(f"Reason: Unphysical. State is outside the Bloch sphere at t={t}, |r|² ≈ {r_sq_t:.3f} > 1.")
            return False

    print("Result: This evolution appears to be physically valid.")
    return True

def main():
    """
    Main function to analyze all plots and find the valid one.
    """
    # Data visually estimated from the plots
    all_data = {
        'A': {'t0': {'sz': 0.4, 'sp_mag': 0.7, 'S': 0.0}, 't_other': [(2.0, {'sz': 0.0, 'sp_mag': 0.85})]},
        'B': {'t0': {'sz': 0.5, 'sp_mag': 0.7, 'S': 0.0}},
        'C': {'invalid_values': [
            {'t': 1.5, 'name': 'sz', 'value': 1.6},
            {'t': 3.0, 'name': 'S', 'value': -1.2}
        ]},
        'D': {'t0': {'sz': 0.5, 'sp_mag': 0.7, 'S': 0.0}},
        'E': {'t0': {'sz': 0.5, 'sp_mag': 0.7, 'S': 0.0}},
        'F': {'t0': {'sz': 0.6, 'sp_mag': 0.4, 'S': 0.0}, 't_other': [
              (0.5, {'sz': 0.5, 'sp_mag': 0.42}), 
              (1.0, {'sz': 0.7, 'sp_mag': 0.3})]
        },
    }

    valid_plot = None
    for plot_name in sorted(all_data.keys()):
        is_valid = check_plot(plot_name, all_data[plot_name])
        if is_valid:
            valid_plot = plot_name
        print("\n")
    
    print("="*30)
    print(f"CONCLUSION: The only physically valid quantum evolution is shown in diagram {valid_plot}.")
    print("="*30)

if __name__ == "__main__":
    main()