def solve_einstein_citation():
    # Step 1: Identify the scientist whose lectures were cited.
    # In his 1905 paper "Über die von der molekularkinetischen Theorie der Wärme geforderte Bewegung
    # von in ruhenden Flüssigkeiten suspendierten Teilchen", Einstein's work on Brownian motion
    # was heavily based on the principles of statistical mechanics, pioneered by Ludwig Boltzmann.
    # Einstein's work can be seen as a direct application and validation of Boltzmann's statistical theories.
    cited_lecturer = "Ludwig Boltzmann"

    # Step 2: Print the answer to the user's question.
    print(f"In his 1905 paper on Brownian motion, Albert Einstein's work was fundamentally based on the principles from the lectures and theories of: {cited_lecturer}.")
    print("-" * 20)

    # Step 3: Display the key equation from the paper.
    # One of the paper's core results is the formula for the diffusion coefficient (D).
    # This equation provided a theoretical way to calculate Avogadro's number (N_A) by observing Brownian motion.
    print("A key equation from the paper is the formula for the diffusion coefficient (D):")
    print("D = (R * T) / (N_A * 6 * \u03c0 * \u03b7 * a)")
    print("Where:")
    print("  R = Gas constant")
    print("  T = Absolute temperature")
    print("  N_A = Avogadro's number")
    print("  \u03b7 = Viscosity of the fluid")
    print("  a = Radius of the particle")
    print("-" * 20)
    
    # Step 4: Output the numbers from the final equation as requested.
    number_in_equation = 6
    print("The numerical constant explicitly mentioned in the denominator of this equation is:")
    print(number_in_equation)

solve_einstein_citation()