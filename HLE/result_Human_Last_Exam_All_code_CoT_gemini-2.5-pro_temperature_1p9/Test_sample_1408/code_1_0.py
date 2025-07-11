import math

def solve():
    """
    Calculates the maximum overhang for three stacked cubes, considering rotations.
    """
    
    # List of all 8 rotation configurations for (Block1, Block2, Block3)
    # 'N' for Normal, 'R' for Rotated
    configs = [
        'NNN', 'NNR', 'NRN', 'NRR',
        'RNN', 'RNR', 'RRN', 'RRR'
    ]

    max_overhang = 0.0
    best_config_details = {}

    sqrt2 = math.sqrt(2)

    for config in configs:
        b1_rot, b2_rot, b3_rot = (c == 'R' for c in config)

        # Determine displacements d1, d2 based on rotation of supporting blocks
        # d1 = c1 - c2, where c1 is CM of B1, c2 is CM of B2
        d1 = 1/sqrt2 if b2_rot else 1/2
        # d2 = (c1+c2)/2 - c3, where c3 is CM of B3
        d2 = 1/sqrt2 if b3_rot else 1/2

        # Calculate CM positions based on the equilibrium conditions
        c1 = d1/2 + d2/3
        c2 = c1 - d1
        c3 = (c1 + c2)/2 - d2

        # Determine half-widths for each block
        w1 = 1/sqrt2 if b1_rot else 1/2
        w2 = 1/sqrt2 if b2_rot else 1/2
        w3 = 1/sqrt2 if b3_rot else 1/2

        # Calculate the overhang for the current configuration
        overhang = max(c1 + w1, c2 + w2, c3 + w3)

        if overhang > max_overhang:
            max_overhang = overhang
            
            # For expressing the overhang in the form (A + B*sqrt(2))/D
            # We work with integer arithmetic on 1 and sqrt(2) components as much as possible
            # to maintain precision. Common denominator for c_i will be 12.
            # d1_num is the numerator of d1 in terms of 1/sqrt(2) or 1
            if b2_rot: d1_val = (0, 1) # 1/sqrt(2) -> 0*1 + 1*1/sqrt(2)
            else:      d1_val = (1, 0) # 1/2
            
            if b3_rot: d2_val = (0, 1)
            else:      d2_val = (1, 0)
            
            # c1 = d1/2 + d2/3. In terms of (coeff of 1, coeff of 1/sqrt(2))
            c1_val = (d1_val[0]/2 + d2_val[0]/3, d1_val[1]/2 + d2_val[1]/3)
            # overhang is c1+w1 in all winning cases. Let's check max term.
            
            # Find which block's edge determines the overhang
            overhangs = {'B1': c1 + w1, 'B2': c2 + w2, 'B3': c3 + w3}
            max_block = max(overhangs, key=overhangs.get)
            
            # Denominator is 12, as it's the LCM of denominators in c_i calculations
            D = 12
            
            if max_block == 'B1':
                # Overhang = c1 + w1 = d1/2 + d2/3 + w1
                # = (1/2 if not b2_rot else 1/(2*sqrt2)) + (1/3 if not b3_rot else 1/(3*sqrt2)) + (1/2 if not b1_rot else 1/sqrt2)
                if not b1_rot and not b2_rot and not b3_rot: #NNN
                     A, B = 11, 0
                elif b1_rot and not b2_rot and not b3_rot: #RNN
                     A, B = 5, 6
                elif not b1_rot and b2_rot and not b3_rot: #NRN
                     A, B = 8, 3
                elif not b1_rot and not b2_rot and b3_rot: #NNR
                     A, B = 9, 2
                elif b1_rot and b2_rot and not b3_rot: #RRN
                     A, B = 2, 9
                elif b1_rot and not b2_rot and b3_rot: #RNR
                     A, B = 3, 8
                elif not b1_rot and b2_rot and b3_rot: #NRR
                     A, B = 6, 5
                elif b1_rot and b2_rot and b3_rot: #RRR
                     A, B = 0, 11
            # In all maximal overhang cases, the top block determines the overhang.
            # No need to code c2+w2, c3+w3 here.

            # Convert to final (a, b, c) format
            # overhang = (A + B*sqrt(2))/12 = (a + sqrt(b))/(1+c)
            #
            # The denominator D=12 implies 1+c=12, so c=11.
            # The numerator is A + B*sqrt(2) = A + sqrt(2*B^2)
            # So a=A, b=2*B^2. If B=0, b=0. If A=0, a=0.
            
            final_a = A
            final_b = 2 * B**2 if B != 0 else 0
            final_c = D - 1
            
            best_config_details = {
                "config": config,
                "a": final_a,
                "b": final_b,
                "c": final_c,
            }

    # Output the final answer
    a = best_config_details["a"]
    b = best_config_details["b"]
    c = best_config_details["c"]
    
    print(f"The maximal overhang is achieved with the configuration '{best_config_details['config']}'.")
    print(f"The overhang value is ({a} + \u221A{b}) / (1 + {c}).")
    print("The requested integers a, b, c are:")
    print(f"{a} {b} {c}")

solve()