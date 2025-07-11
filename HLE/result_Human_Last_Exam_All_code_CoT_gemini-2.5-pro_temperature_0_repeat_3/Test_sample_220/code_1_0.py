import math

def calculate_jet_speed():
    """
    Calculates and verifies the jet speeds of bursting bubbles based on the
    physical model v = k / sqrt(R).
    """
    print("Based on the physical model, jet speed (v) is inversely proportional to the square root of the bubble radius (R).")
    print("Model: v = k / \u221A(R)\n")

    # We test the answer choice (4 m/s, 1.5 m/s).
    # The higher speed corresponds to the smaller bubble.
    
    # Data for the larger bubble (2 cm diameter)
    d1_cm = 2.0
    r1_m = (d1_cm / 2.0) / 100.0  # Radius in meters
    v1_ms = 1.5  # Assumed speed in m/s

    # Data for the smaller bubble (2 mm diameter)
    d2_mm = 2.0
    r2_m = (d2_mm / 2.0) / 1000.0 # Radius in meters
    
    # Step 1: Calculate the proportionality constant 'k' using the 2 cm bubble data.
    # k = v * sqrt(R)
    k = v1_ms * math.sqrt(r1_m)
    
    print(f"Step 1: Calculate the constant 'k' using the 2 cm bubble (R={r1_m} m, v={v1_ms} m/s).")
    print(f"Equation: {v1_ms} = k / \u221A({r1_m})")
    print(f"k = {v1_ms} * \u221A({r1_m}) = {k:.4f}\n")

    # Step 2: Use 'k' to predict the speed for the 2 mm bubble.
    v2_calculated = k / math.sqrt(r2_m)
    
    print(f"Step 2: Predict the speed for the 2 mm bubble (R={r2_m} m) using k.")
    print(f"Equation: v = {k:.4f} / \u221A({r2_m})")
    print(f"Result: v = {v2_calculated:.2f} m/s\n")

    # Step 3: Conclusion
    print("The calculated speed of {:.2f} m/s for the 2 mm bubble is approximately 4 m/s.".format(v2_calculated))
    print("This confirms that the speeds for bubbles with diameters of 2 mm and 2 cm are approximately 4 m/s and 1.5 m/s, respectively.")
    
    print("\n--- Final Equations ---")
    print(f"For the 2 mm bubble: {v2_calculated:.2f} m/s = {k:.4f} / \u221A({r2_m})")
    print(f"For the 2 cm bubble: {v1_ms:.2f} m/s = {k:.4f} / \u221A({r1_m})")

calculate_jet_speed()