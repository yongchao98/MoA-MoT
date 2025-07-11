def solve_chemistry_problem():
    """
    This script simulates the chemical reactions described to determine the final barium salt.
    """

    print("Step 1: Mixing aqueous solutions of barium chloride (BaCl2) and silver nitrate (AgNO3).")
    print("A double displacement reaction occurs:")
    print("1 BaCl2(aq) + 2 AgNO3(aq) -> 1 Ba(NO3)2(aq) + 2 AgCl(s)")
    print("Based on solubility rules, Silver Chloride (AgCl) is an insoluble solid and precipitates.")
    print("Barium Nitrate (Ba(NO3)2) is soluble and remains dissolved in the water.")
    print("--> Barium salt in the flask is: Barium Nitrate (Ba(NO3)2)\n")

    print("Step 2: The mixture is dried using freeze-drying.")
    print("The water is removed, leaving behind the solid products.")
    print("--> Flask now contains solid Barium Nitrate (Ba(NO3)2) and solid Silver Chloride (AgCl).\n")

    print("Step 3: Ammonia (NH3) is added and mixed.")
    print("Ammonia reacts with solid Silver Chloride to form a soluble complex, diamminesilver(I) chloride.")
    print("Reaction: AgCl(s) + 2 NH3(aq) -> [Ag(NH3)2]+(aq) + Cl-(aq)")
    print("Barium Nitrate does not react with ammonia.")
    print("--> The Barium salt in the flask is still solid Barium Nitrate (Ba(NO3)2).\n")

    print("Step 4: The ammonia is evaporated using freeze-drying.")
    print("Removing the ammonia reverses the complex formation from Step 3.")
    print("The diamminesilver(I) complex breaks down, and Silver Chloride (AgCl) precipitates as a solid again.")
    print("The Barium Nitrate is unaffected by this step.")
    print("--> The final contents are solid Barium Nitrate (Ba(NO3)2) and solid Silver Chloride (AgCl).\n")

    print("Conclusion:")
    final_barium_salt = "Barium Nitrate"
    final_formula = "Ba(NO3)2"
    print(f"The barium salt that is in the flask after all reactions is {final_barium_salt}.")
    print(f"Its chemical formula is: {final_formula}")

solve_chemistry_problem()