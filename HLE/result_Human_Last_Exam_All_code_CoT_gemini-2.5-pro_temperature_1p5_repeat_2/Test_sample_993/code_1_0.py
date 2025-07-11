def solve_paleoclimate_question():
    """
    Explains the reasoning behind the isotopic composition of sulfate
    in lakebeds and determines the correct answer choice.
    """

    explanation = """
    Step-by-step analysis to determine the factors responsible for low levels of the 18O isotope in sulfate deposits:

    1.  **Isotope Source:** The oxygen in lake sulfate (SO4) comes from the lake water (H2O). Therefore, the isotopic composition of the lake water (the ratio of heavy 18O to light 16O) is the primary driver.

    2.  **Evaporation's Effect:** Evaporation preferentially removes the lighter 16O isotope, leaving the remaining lake water enriched in the heavier 18O. Dry and warm climates maximize evaporation, leading to HIGH 18O levels in the lake water and the resulting sulfate.

    3.  **Precipitation's Effect:** Precipitation (rain/snow) is naturally depleted in 18O. A wet climate means the lake is dominated by inflow from precipitation, which keeps the 18O levels LOW. Furthermore, cold climates cause precipitation to be even more depleted in 18O.

    4.  **Lake Level as an Indicator:** High lake levels indicate that water inflow from precipitation is greater than water loss from evaporation. This is characteristic of a wet climate and creates a large reservoir that is less affected by evaporation, helping to keep 18O levels low.

    5.  **Conclusion:** To achieve LOW 18O levels in sulfate deposits, we need conditions that minimize evaporation and maximize inflow of 18O-depleted water. This corresponds to a WET, COLD climate, which would naturally lead to HIGH lake levels.
    """

    answer_choice = 'G'
    answer_text = 'Wet, cold climate with high lake levels'

    print(explanation)
    print(f"The correct combination of factors is: {answer_choice}. {answer_text}")

solve_paleoclimate_question()