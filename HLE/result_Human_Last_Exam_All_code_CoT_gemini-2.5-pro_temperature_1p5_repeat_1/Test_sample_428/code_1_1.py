import sympy

def solve_particle_annihilation():
    """
    Calculates the expected time until only one particle remains.
    The user will be prompted to enter the positive integer values for
    N1, M1, N2, and M2.
    """
    try:
        N1_in = input("Enter the value for N1: ")
        M1_in = input("Enter the value for M1: ")
        N2_in = input("Enter the value for N2: ")
        M2_in = input("Enter the value for M2: ")

        N1 = sympy.Symbol('N1')
        M1 = sympy.Symbol('M1')
        N2 = sympy.Symbol('N2')
        M2 = sympy.Symbol('M2')
        
        # Initial positions of the 5 particles
        x = [0,
             N1,
             N1 + M1,
             N1 + M1 + N2,
             N1 + M1 + N2 + M2]

        # 1. Calculate the sum of squared differences S(0)
        s0 = 0
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                s0 += (x[j] - x[i])**2
        
        s0 = sympy.simplify(s0)
        
        # 2. Calculate the expected time of the first collision, E[tau1]
        # Using the formula E[tau1] = (d1+d2)(d3+d4)/2
        # where d1=N1, d2=M1, d3=N2, d_4=M2
        e_tau1 = (N1 + M1) * (N2 + M2) / 2
        
        # 3. Calculate the expected time of the second collision, E[tau2]
        # From the martingale equation: 8 * E[tau1] + 12 * E[tau2] = S(0)
        e_tau2 = (s0 - 8 * e_tau1) / 12
        e_tau2_simplified = sympy.simplify(e_tau2)

        # Substitute user provided values
        N1_val = int(N1_in)
        M1_val = int(M1_in)
        N2_val = int(N2_in)
        M2_val = int(M2_in)
        
        s0_val = s0.subs({N1: N1_val, M1: M1_val, N2: N2_val, M2: M2_val})
        e_tau1_val = e_tau1.subs({N1: N1_val, M1: M1_val, N2: N2_val, M2: M2_val})
        e_tau2_val = e_tau2_simplified.subs({N1: N1_val, M1: M1_val, N2: N2_val, M2: M2_val})
        
        # Output the components of the final calculation
        # The final equation is E[tau] = (S(0) - 8 * E[tau1]) / 12
        print(f"The calculation is based on the equation: E[tau] = (S(0) - 8 * E[tau1]) / 12")
        print(f"Given N1={N1_val}, M1={M1_val}, N2={N2_val}, M2={M2_val}:")
        print(f"Initial sum of squared differences S(0) = {s0_val}")
        print(f"Expected time of the first collision E[tau1] = {e_tau1_val}")
        print(f"Therefore, the expected time of the second collision E[tau] = ({s0_val} - 8 * {e_tau1_val}) / 12")
        print(f"E[tau] = {e_tau2_val}")

    except (ValueError, TypeError) as e:
        print(f"Invalid input. Please enter positive integers. Error: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_particle_annihilation()