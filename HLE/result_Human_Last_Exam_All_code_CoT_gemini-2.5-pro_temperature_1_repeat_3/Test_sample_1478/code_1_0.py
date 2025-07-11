#
# A Python program to simulate the Titan 6-bit architecture
# for calculating gravitational force, demonstrating its feasibility.
#

def run_titan_simulation():
    """
    This function simulates the calculation on the Titan architecture
    and prints the step-by-step process.
    """
    BIT_LIMIT = 6
    MAX_INT = 2**BIT_LIMIT - 1  # Max value for a 6-bit integer is 63

    print("--- Titan Feasibility Study for Gravity Calculation ---")
    print(f"Constraint: All numerators/denominators must be <= {MAX_INT}\n")
    print("Objective: Calculate Mass M = ρ * (4/3) * π * r³")
    print("The fractional part of this calculation is: 12 * (4/3) * (22/7) * 8\n")

    try:
        # Step 1: Start with the density factor, 12
        # MOV AX, 12
        ax_val = 12
        print(f"1. MOV AX, {ax_val}")
        print(f"   > AX now holds the value {ax_val}.")

        # Step 2: Multiply by (4/3)
        # MUL AX, 4/3
        mul_num, mul_den = 4, 3
        print(f"\n2. MUL AX, {mul_num}/{mul_den}")
        
        res_num = ax_val * mul_num
        res_den = mul_den
        print(f"   > Calculating ({ax_val} * {mul_num}) / {res_den} = {res_num}/{res_den}")
        
        # A REDUCE instruction would simplify the fraction
        simplified_num = res_num // res_den
        simplified_den = 1
        print(f"   > REDUCE AX. Result is {simplified_num}/{simplified_den}.")
        
        if simplified_num > MAX_INT or simplified_den > MAX_INT:
            raise ValueError("Value exceeds 6-bit limit.")
        ax_val = simplified_num
        print(f"   > OK. AX now holds {ax_val}.")

        # Step 3: Multiply by Pi (approximated as 22/7)
        # MUL AX, 22/7
        mul_num, mul_den = 22, 7
        print(f"\n3. MUL AX, {mul_num}/{mul_den}")
        
        res_num = ax_val * mul_num
        print(f"   > Calculating ({ax_val} * {mul_num}) / {mul_den} = {res_num}/{mul_den}")

        if res_num > MAX_INT:
            print(f"   > FAILURE: Numerator {res_num} is greater than {MAX_INT}.")
            print("   > Attempting expansion trick as per instructions: 22/7 = 3 + 1/7.")
            print(f"   > Expression becomes: {ax_val} * (3 + 1/7) = ({ax_val} * 3) + ({ax_val}/7)")
            
            part1_num = ax_val * 3
            part2_num = ax_val
            part2_den = 7
            print(f"   > Part 1: {ax_val}*3 = {part1_num}. Part 2: {part2_num}/{part2_den}.")
            print(f"   > The expression held in AX is now: ({part1_num}/1 + {part2_num}/{part2_den})")
            
            # Step 4: This new expression must be multiplied by the next term (r³ factor, 8)
            print("\n4. Next, multiply the expression in AX by the r³ factor, 8.")
            print("   > This requires multiplying each term in the expression by 8:")
            print(f"   > ({part1_num} * 8) + (({part2_num}/{part2_den}) * 8)")
            
            final_res_num = part1_num * 8
            print(f"   > Calculating first term: {part1_num} * 8 = {final_res_num}")
            
            if final_res_num > MAX_INT:
                print(f"   > FAILURE: Resulting numerator {final_res_num} is greater than {MAX_INT}.")
                raise ValueError("Calculation failed due to 6-bit integer overflow.")
            else:
                # This part would not be reached, but included for completeness
                print("   > Calculation succeeded.")

    except ValueError as e:
        print(f"\n--- CALCULATION HALTED ---")
        print(f"Reason: {e}")
        print("\nConclusion: The calculation is not feasible on the Titan architecture.")
        print("The magnitude of intermediate products (e.g., 16 * 22 = 352) exceeds the")
        print("6-bit integer limit of 63. Furthermore, the expansion-and-distribution strategy")
        print("also fails because the resulting components are still too large to be multiplied")
        print("in subsequent steps (e.g., 48 * 8 = 384). The architecture's constraints")
        print("are too restrictive for the numbers involved in this physics problem.")
        print("\nFinal Answer is therefore N0, as the calculation cannot be completed.")
        print("<<<N0>>>")

if __name__ == '__main__':
    run_titan_simulation()