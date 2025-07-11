import sys

def solve_maxwell_demon():
    """
    This function explains the physics behind the Maxwell's demon apparatus and determines the critical parameter.
    """
    # The problem describes a physical implementation of Maxwell's demon, where a one-way door is intended to separate gas molecules.
    # The goal is to trap all gas molecules on one side, which seemingly violates the Second Law of Thermodynamics.

    # According to physicist Richard Feynman's analysis of such a system (a Brownian ratchet), the door itself must be considered.
    # If the door is at the same temperature as the gas, it will have its own thermal energy and will jiggle randomly.
    # This random motion will cause the door to occasionally open and let molecules pass in the reverse direction, preventing a perfect separation.

    # For the door to act as a perfect one-way gate, it must have no thermal fluctuations.
    # The physical parameter that measures thermal energy is Temperature.
    # The state of zero thermal energy is reached at absolute zero.

    # Therefore, we can write a conceptual equation for the ideal condition of the door:
    # T_door = 0 Kelvin
    
    equation_text = "T_door = 0"
    
    print("For the one-way door to perfectly trap all gas, it must not be affected by the random thermal motion of the system.")
    print("This requires the door itself to have no thermal energy.")
    print("The physical parameter that corresponds to thermal energy is Temperature.")
    print("The condition for the door to work perfectly can be expressed with the following conceptual equation, where T_door is the door's temperature in Kelvin:")
    
    # Extracting and printing the numbers from the equation as requested
    numbers_in_equation = []
    current_number = ""
    for char in equation_text:
        if char.isdigit() or char == '.':
            current_number += char
        elif current_number:
            numbers_in_equation.append(current_number)
            current_number = ""
    if current_number:
        numbers_in_equation.append(current_number)

    print("T_door = ", end="")
    for number in numbers_in_equation:
        print(f"{number}", end=" ")
    print("Kelvin")

    print("\nThus, Temperature is the critical experimental parameter that needs to be controlled.")

solve_maxwell_demon()