import math

def solve_shannon_capacity():
    """
    Calculates the Shannon capacity of G⊠H based on the derivation.
    G is K_m with a C5 removed.
    H is K_n with a C4 removed.
    """
    print("This program calculates the Shannon capacity of G⊠H.")
    print("Please provide the natural numbers n and m.")

    while True:
        try:
            m_str = input("Enter the number of vertices m (must be >= 5): ")
            m = int(m_str)
            if m < 5:
                print("Error: m must be greater than or equal to 5.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter an integer.")

    while True:
        try:
            n_str = input("Enter the number of vertices n (must be >= 4): ")
            n = int(n_str)
            if n < 4:
                print("Error: n must be greater than or equal to 4.")
                continue
            break
        except ValueError:
            print("Invalid input. Please enter an integer.")

    print("\n--- Derivation ---")

    # Step 1: Analyze H
    print("1. Capacity of H:")
    print("   H is K_n with a C_4 removed. Its complement H_bar is a C_4 and n-4 isolated vertices.")
    print("   H is a perfect graph, so its capacity c(H) is its independence number alpha(H).")
    print("   alpha(H) = omega(H_bar) = 2.")
    c_H = 2
    print(f"   Therefore, c(H) = {c_H}\n")

    # Step 2: Analyze G
    print("2. Capacity of G:")
    print("   G is K_m with a C_5 removed. G is not perfect.")
    print("   G can be constructed as the join of a C_5 and K_{m-5}.")
    print("   The capacity is c(G) = Theta(C_5) + m - 5.")
    sqrt_5 = math.sqrt(5)
    c_G = sqrt_5 + m - 5
    print(f"   c(G) = sqrt(5) + {m} - 5 = {sqrt_5:.4f} + {m-5} = {c_G:.4f}\n")

    # Step 3: Combine for G⊠H
    print("3. Capacity of the Strong Product G⊠H:")
    print("   Since H is a perfect graph, the capacity of the strong product is multiplicative:")
    print("   c(G⊠H) = c(G) * c(H)")

    final_capacity = c_G * c_H

    # Final result and equation
    print("\n--- Final Equation and Result ---")
    print(f"The Shannon capacity is given by the equation:")
    print(f"Capacity = c(G) * c(H)")
    print(f"Capacity = (sqrt(5) + {m} - 5) * {c_H}")
    print(f"The numbers in the final equation are 2, 5, and 5, with m={m}.")
    print(f"Result for m={m} and n={n}:")
    print(f"Capacity = ({sqrt_5:.4f} + {m-5}) * {c_H} = {final_capacity:.4f}")

if __name__ == '__main__':
    solve_shannon_capacity()
