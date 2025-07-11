import sys

def calculate_ni_ce_ratio():
    """
    Calculates and explains the ideal Ni/Ce molar ratio based on findings
    from scientific literature for Ni-Ceria catalysts.
    """
    # Based on experimental literature, a common optimal loading is 10 mol% Ni.
    # This means in a sample of 100 total moles of metal (Ni + Ce), 10 are Ni.
    optimal_ni_moles = 10
    total_moles = 100

    # The rest of the moles are Cerium.
    ce_moles = total_moles - optimal_ni_moles

    # The ratio is Ni moles divided by Ce moles.
    ni_ce_ratio = optimal_ni_moles / ce_moles

    print("The ideal Ni/Ce ratio is determined experimentally. A widely cited optimal value is around 10 mol% Ni.")
    print("This balances the number of active sites with nanoparticle dispersion and strong metal-support interaction.")
    print("\nHere is the calculation to find the direct Ni to Ce molar ratio from this value:")
    print(f"1. Assume a total of {total_moles} moles of metal atoms (Ni + Ce).")
    print(f"2. A 10 mol% Ni content means there are {optimal_ni_moles} moles of Ni.")
    print(f"3. The number of Ce moles is therefore Total - Ni = {total_moles} - {optimal_ni_moles} = {ce_moles} moles.")
    print(f"\n4. The final equation for the Ni/Ce ratio is:")
    print(f"   Ni moles   /   Ce moles   =   Ratio")
    print(f"   {optimal_ni_moles:<10} /   {ce_moles:<10} =   {ni_ce_ratio:.3f}")
    print(f"\nTherefore, the ideal Ni/Ce molar ratio is approximately {ni_ce_ratio:.3f}.")

if __name__ == '__main__':
    calculate_ni_ce_ratio()
