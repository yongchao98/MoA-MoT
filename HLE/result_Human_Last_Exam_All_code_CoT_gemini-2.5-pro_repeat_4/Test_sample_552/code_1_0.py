def find_cited_lecturer():
    """
    This function identifies and returns the name of the physicist whose lectures
    on statistical mechanics and gas theory were a primary influence on Albert
    Einstein's 1905 paper on Brownian motion.

    While the paper itself lacks a formal bibliography, historical analysis confirms
    that Einstein's work is built upon the foundations laid by this scientist.
    """
    # The physicist is Ludwig Boltzmann. His work on the molecular-kinetic theory
    # of heat was the theoretical basis for Einstein's paper.
    influential_physicist = "Ludwig Boltzmann"
    return influential_physicist

def main():
    """
    Main function to execute the logic and print the result.
    """
    lecturer_name = find_cited_lecturer()
    print("In his 1905 paper on Brownian motion, Albert Einstein built upon the theoretical framework developed by a specific physicist.")
    print("The lectures he studied, which were foundational to his paper, were by:")
    print(f"\n{lecturer_name}")

if __name__ == "__main__":
    main()