def design_liquid_crystal():
    """
    This script models the design of a liquid crystal molecule to achieve a
    transition temperature near room temperature, based on a given set of rules.
    """
    # The general structure is C(n)H(2n+1)-Ph-CN
    # Target transition temperature is near room temperature (e.g., 20-25 degrees Celsius)

    print("Step 1: Define the target and initial parameters.")
    target_temp_min = 20.0
    target_temp_max = 25.0
    # Rule F: Start with a pentyl chain (n=5)
    n = 5
    print(f"Target temperature range: {target_temp_min}°C - {target_temp_max}°C")
    print(f"Starting with initial chain length n = {n}\n")

    # Heuristic formula based on rule F: "If transition temperature is too high, increase chain length"
    # This model uses the formula: T(n) = 35 - 5 * (n - 5)
    def calculate_transition_temp(n_val):
        return 35.0 - 5.0 * (n_val - 5.0)

    print("Step 2: Iteratively search for the optimal chain length 'n'.")
    while True:
        # Calculate temperature for the current 'n'
        current_temp = calculate_transition_temp(n)

        print(f"Testing n = {n}...")
        print(f"  - Molecular structure: C({n})H({2*n+1})-Ph-CN")
        print(f"  - Estimated transition temperature: {current_temp:.1f}°C")

        # Check if temperature is in the target range
        if target_temp_min <= current_temp <= target_temp_max:
            print("\nSuccess! Temperature is within the target range.")
            break

        # Rule F: Tuning logic
        if current_temp > target_temp_max:
            print("  - Temperature is too high. Increasing chain length 'n'.\n")
            n += 1
        elif current_temp < target_temp_min:
            print("  - Temperature is too low. Decreasing chain length 'n'.\n")
            n -= 1
        
        # Safety break to prevent infinite loops
        if n > 20 or n < 1:
            print("\nCould not find a suitable chain length within the range n=1 to 20.")
            return

    print("\nStep 3: Final proposed molecular design.")
    final_n = n
    final_h = 2 * final_n + 1
    final_temp = current_temp

    print("The final proposed molecule has the following properties:")
    print(f" - Optimal chain length n: {final_n}")
    print(f" - Estimated transition temperature: {final_temp:.1f}°C")
    
    # Output each number in the final equation
    print("\nFinal Equation (Molecular Formula):")
    print(f"C({final_n})H({final_h})-Ph-CN")


# Execute the design process
design_liquid_crystal()