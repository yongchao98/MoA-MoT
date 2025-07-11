import re

def solve_chemistry_problem():
    """
    This function simulates the chemical reactions described to find the final barium salt.
    """
    # Step 1: Reaction of Barium Chloride and Silver Nitrate
    # BaCl2(aq) + 2AgNO3(aq) -> Ba(NO3)2(aq) + 2AgCl(s)
    barium_salt = {"name": "Barium Nitrate", "formula": "Ba(NO3)2"}
    print(f"Step 1: Mixing Barium Chloride and Silver Nitrate yields {barium_salt['name']} and Silver Chloride.")

    # Step 2 & 3: Drying and adding Ammonia
    # Barium Nitrate does not react with Ammonia.
    print(f"Step 2: After drying and adding ammonia, the barium salt remains {barium_salt['name']}.")

    # Step 4: Final Drying
    # The Barium Nitrate is unaffected by the removal of ammonia.
    print(f"Step 3: After the final drying step, the barium salt is still {barium_salt['name']}.")

    print("-" * 30)
    print(f"The final barium salt in the flask is {barium_salt['name']}.")
    
    formula = barium_salt['formula']
    print(f"The chemical formula for the final barium salt is: {formula}")
    
    # Extracting numbers from the formula
    numbers = re.findall(r'\d+', formula)
    
    print("The numbers in the final equation's formula are:")
    for num in numbers:
        print(num)

solve_chemistry_problem()