import math

def is_prime(n):
    """Checks if a number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def solve():
    """
    Solves the circuit problem by searching for resistor values that satisfy
    all constraints and calculating the maximum possible current through R3.
    """
    R2 = 6
    max_current = 0.0
    best_config = None

    # Search for valid (R1, R3) pairs.
    # Analysis shows R1 must be > 8, so we can start the search there.
    # A reasonable upper bound is chosen for the search space.
    for R1 in range(9, 50):
        # Per R3 - R1 > 2, R3 must be at least R1 + 3.
        for R3_candidate in range(R1 + 3, 100):
            if not is_prime(R3_candidate):
                continue
            
            R3 = R3_candidate
            
            # Check the z(C)=6 condition. 
            # This is equivalent to R2=6 being the value that is "farthest" from the mean of the other two.
            d_for_R1 = abs(R1 - (R2 + R3) / 2.0)
            d_for_R2 = abs(R2 - (R1 + R3) / 2.0)
            d_for_R3 = abs(R3 - (R1 + R2) / 2.0)
            
            if d_for_R2 > d_for_R1 and d_for_R2 > d_for_R3:
                # This is a valid (R1, R3) pair. Now calculate currents for both scenarios.
                
                # --- Scenario A: R1, R2, R3 in parallel ---
                # When R2 fails, circuit is R1 || R3. Voltage across them is V_prime = 26V.
                # Total source current I = I_prime_1 + I_prime_3 = 26/R1 + 26/R3.
                # With R2 intact, this source current I is split across R1, R2, R3.
                # The equivalent resistance is R_eq = 1 / (1/R1 + 1/R2 + 1/R3).
                # The voltage is V_intact = I * R_eq.
                # The current through R3 is I3 = V_intact / R3.
                # A simplified formula is used below.
                
                num_A = 26.0 * (R1 + R3) * R2
                den_A = float(R1 * R2 + R1 * R3 + R2 * R3)
                current_A = num_A / den_A

                if current_A > max_current:
                    max_current = current_A
                    best_config = {
                        "scenario": "A", "R1": R1, "R3": R3, 
                        "current": current_A, "num": num_A, "den": den_A
                    }
                
                # --- Scenario B: R1 in series with (R2 || R3) ---
                # When R2 fails, circuit is R1 in series with R3.
                # V_prime_3 = I_source * R3 = 26 => I_source = 26 / R3.
                # With R2 intact, source current I_source splits between R2 and R3.
                # Current through R3 is I3 = I_source * (R2 / (R2 + R3)).
                num_B = 26.0 * R2
                den_B = float(R3 * (R2 + R3))
                current_B = num_B / den_B

                if current_B > max_current:
                    max_current = current_B
                    best_config = {
                        "scenario": "B", "R1": R1, "R3": R3,
                        "current": current_B, "num": num_B, "den": den_B
                    }

    # Print the final results based on the best configuration found
    if best_config:
        R1 = best_config["R1"]
        R3 = best_config["R3"]
        scenario = best_config["scenario"]
        
        print(f"The maximum possible current is achieved with R1 = {R1} ohms and R3 = {R3} ohms.")
        
        if scenario == "A":
            print("This occurs in Scenario A, where all three resistors are connected in parallel.")
            print("\nThe current through R3 (I_3) in the intact circuit is calculated as:")
            print("I_3 = (26 * (R_1 + R_3) * R_2) / (R_1*R_2 + R_1*R_3 + R_2*R_3)")
            print("\nPlugging in the numerical values:")
            print(f"I_3 = (26 * ({R1} + {R3}) * {R2}) / ({R1}*{R2} + {R1}*{R3} + {R2}*{R3})")
            num = int(best_config["num"])
            den = int(best_config["den"])
            print(f"I_3 = {num} / {den}")
            print(f"Maximum possible current = {best_config['current']} A")
        else: # scenario B
            print("This occurs in Scenario B, where R1 is in series with the parallel R2 and R3.")
            print("\nThe current through R3 (I_3) in the intact circuit is calculated as:")
            print("I_3 = (26 * R_2) / (R_3 * (R_2 + R_3))")
            print("\nPlugging in the numerical values:")
            print(f"I_3 = (26 * {R2}) / ({R3} * ({R2} + {R3}))")
            num = int(best_config["num"])
            den = int(best_config["den"])
            print(f"I_3 = {num} / {den}")
            print(f"Maximum possible current = {best_config['current']} A")

solve()
<<<13.388059701492538>>>