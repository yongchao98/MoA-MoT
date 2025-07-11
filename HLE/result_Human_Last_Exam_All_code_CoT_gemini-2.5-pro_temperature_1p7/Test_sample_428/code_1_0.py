import sys

# Define the positive integer variables N1, M1, N2, M2.
# You can change these placeholder values to your specific numbers.
try:
    N1 = 1
    M1 = 2
    N2 = 3
    M2 = 4

    # Check if variables were passed as command-line arguments
    if len(sys.argv) == 5:
        N1 = int(sys.argv[1])
        M1 = int(sys.argv[2])
        N2 = int(sys.argv[3])
        M2 = int(sys.argv[4])
        
    if N1 <= 0 or M1 <= 0 or N2 <= 0 or M2 <= 0:
        raise ValueError("The values for N1, M1, N2, M2 must be positive integers.")

except (ValueError, IndexError):
    print("Using placeholder values. To use your own values, either edit the script or run as:")
    print("python your_script_name.py N1 M1 N2 M2")
    print("For example: python your_script_name.py 1 2 3 4")
    N1 = 1
    M1 = 2
    N2 = 3
    M2 = 4

# The rates of the random walks
lambda1 = 1.0
lambda2 = 2.0

# The formula for the expected time E[tau] until only one particle remains.
# E[tau] = E[tau1] + E[tau2_prime]
# E[tau1] = (N1*M1 + N2*M2 + 0.5*(N1+M2)*(M1+N2) - 0.5*M1*N2) / lambda1
# E[tau2_prime] = (0.5 * M1 * N2) / lambda2
# Adding them up simplifies the expression.

term1 = N1 * M1
term2 = N2 * M2
term3_factor1 = N1 + M2
term3_factor2 = M1 + N2
term3 = 0.5 * term3_factor1 * term3_factor2

expected_tau = term1 + term2 + term3

print("The formula for the expectation of tau is:")
print(f"E[tau] = {N1} * {M1} + {N2} * {M2} + 0.5 * ({N1} + {M2}) * ({M1} + {N2})")
print("\nWith the given values:")
print(f"N1 = {N1}")
print(f"M1 = {M1}")
print(f"N2 = {N2}")
print(f"M2 = {M2}")
print(f"\nThe expectation of tau is: {expected_tau}")
