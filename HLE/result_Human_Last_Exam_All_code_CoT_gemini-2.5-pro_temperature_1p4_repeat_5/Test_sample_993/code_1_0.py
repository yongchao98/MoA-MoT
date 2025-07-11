import textwrap

def solve_paleoclimate_question():
    """
    Analyzes the factors responsible for low levels of the 18O isotope
    in ancient lakebed sulfate deposits and prints the correct answer choice.
    """
    explanation = """
    To determine the factors causing low levels of the 18O isotope in lake sulfate, we need to consider the processes of evaporation and inflow.

    1. Evaporation and Isotopes: The lighter oxygen isotope, 16O, evaporates more easily than the heavier 18O. High rates of evaporation therefore leave the remaining lake water enriched in 18O. To have low levels of 18O, conditions must minimize evaporation. This points to a cold and wet climate.

    2. Inflow and Isotopes: Precipitation that forms in colder climates is naturally depleted in 18O. This is because the heavier 18O rains out more easily as water vapor travels from its ocean source. Therefore, a cold climate provides inflow water that is already low in 18O.

    3. Lake Level: High lake levels indicate a large water volume, which buffers the lake's isotopic composition against the effects of evaporation. A high lake level also implies that inflow is greater than or equal to outflow/evaporation, a condition found in a wet climate.

    Conclusion: The combination of a 'wet' climate (reducing net evaporation), a 'cold' climate (reducing the rate of evaporation and providing 18O-depleted inflow), and 'high lake levels' (buffering against evaporation effects) creates the ideal environment for sulfate deposits with low 18O levels.
    """

    answer_choice = "G"
    full_answer_text = "Wet, cold climate with high lake levels"

    print(textwrap.dedent(explanation).strip())
    print("\n----------------------------------")
    print(f"The correct option is: {answer_choice}. {full_answer_text}")
    print(f"<<<{answer_choice}>>>")

solve_paleoclimate_question()