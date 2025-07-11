import math

def solve_computational_factor():
    """
    This script identifies the original computational factor from the relevant scientific literature.
    """

    print("Step 1: Identify the context and the key scientific papers.")
    print("The problem describes the enthalpy-porosity technique for modeling melting.")
    print(" - The 'prior published simulation-only work' refers to Voller & Prakash (1987).")
    print(" - The 'benchmarked against the melting of...gallium' work refers to Brent, Voller & Reid (1988).\n")

    print("Step 2: Define the Carman-Kozeny source term mentioned in the problem.")
    print("The source term 'S' added to the momentum equation has the form:")
    print("S = -A * u, where A = C * (1-f)^2 / (f^3 + b)\n")
    print("Here, 'C' is the computational factor in question.\n")

    print("Step 3: Find the value of 'C' in each paper.")
    # Value from the prior, simulation-only paper (Voller & Prakash, 1987)
    c_original_base = 1.6
    c_original_exponent = 6
    c_original = c_original_base * (10**c_original_exponent)

    # Value from the later, benchmark paper (Brent, Voller & Reid, 1988)
    c_benchmark_base = 1.6
    c_benchmark_exponent = 3
    c_benchmark = c_benchmark_base * (10**c_benchmark_exponent)

    print(f" - In the prior work (1987), the value used was {c_original_base} x 10^{c_original_exponent}.")
    print(f" - In the later benchmark work (1988), the value was modified to {c_benchmark_base} x 10^{c_benchmark_exponent} to match experiments.\n")

    print("Step 4: Answer the specific question.")
    print("The question asks for the value originally used in the prior published implementation.\n")

    print("Final Answer: The original computational factor 'C' was:")
    # Using 'math.log10' and formatting to output the numbers clearly
    base = c_original / (10**int(math.log10(c_original)))
    exponent = int(math.log10(c_original))
    print(f"C = {base} * 10^{exponent}")

if __name__ == '__main__':
    solve_computational_factor()