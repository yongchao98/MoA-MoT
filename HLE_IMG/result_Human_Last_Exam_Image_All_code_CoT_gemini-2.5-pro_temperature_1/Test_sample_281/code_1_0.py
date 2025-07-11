import sys

def calculate_conductance():
    """
    Calculates the four-terminal conductance G_12,34 for a quantum Hall device
    with a quantum point contact (QPC).
    """
    try:
        m_str = input("Enter the total number of spin-degenerate edge states (M): ")
        M = int(m_str)
        if M <= 0:
            print("Error: M must be a positive integer.")
            return

        n_str = input("Enter the number of reflected edge states by the QPC (N): ")
        N = int(n_str)

        if not (0 <= N <= M):
            print(f"Error: N must be between 0 and M (inclusive). You entered N={N}, M={M}.")
            return

        print("\n--- Calculation ---")
        print("The four-terminal conductance G_12,34 is given by the formula:")
        print("G_12,34 = [M * (M - N) / N] * (e^2/h)\n")

        print(f"For M = {M} and N = {N}:")

        if N == 0:
            # This corresponds to V34 = 0, so the conductance is theoretically infinite.
            print("Since N = 0, the QPC is fully transparent to the inner channels.")
            print("The measured voltage V_34 = V_3 - V_4 would be 0.")
            print("G_12,34 = I / 0 -> infinity")
            print("The conductance is infinite (or a short circuit).")
            # We don't have a numerical answer to return in this case.
            # So no <<<...>>> tag
            return
            
        m_minus_n = M - N
        numerator = M * m_minus_n
        denominator = N
        
        conductance = numerator / denominator

        # Print the equation with substituted values
        print(f"G_12,34 = [{M} * ({M} - {N}) / {N}] * (e^2/h)")
        print(f"G_12,34 = [{M} * {m_minus_n} / {N}] * (e^2/h)")
        print(f"G_12,34 = [{numerator} / {N}] * (e^2/h)")
        print(f"G_12,34 = {conductance:.4f} * (e^2/h)")
        
        # Final answer format
        print(f"\nFinal Answer (in units of e^2/h):")
        print(f"<<<{conductance}>>>")

    except ValueError:
        print("Invalid input. Please enter integers for M and N.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    calculate_conductance()