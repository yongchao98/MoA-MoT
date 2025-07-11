import sys

def solve_paleoclimate_question():
    """
    This script explains the reasoning to determine the climate conditions
    that lead to low levels of the 18O isotope in lakebed sulfate.
    """
    
    print("Step 1: Understand the Proxy")
    print("The oxygen isotopes in sulfate (SO4) deposits reflect the oxygen isotope composition of the lake water (H2O) they formed in.")
    print("Conclusion: Low 18O in sulfate means the ancient lake water had low 18O levels.")
    print("-" * 30)

    print("Step 2: Analyze the Effect of Temperature on Precipitation")
    print("During the water cycle, water with the heavier 18O isotope precipitates more easily than water with the lighter 16O isotope.")
    print("In cold climates, air masses lose a large fraction of their 18O through precipitation before reaching the lake's watershed.")
    print("Conclusion: A 'Cold Climate' is necessary for the inflow (rain/snow) to be depleted in 18O.")
    print("-" * 30)
    
    print("Step 3: Analyze the Effect of Evaporation on the Lake")
    print("Evaporation from a lake surface preferentially removes the lighter 16O isotope.")
    print("This process leaves the remaining lake water enriched in 18O, raising its 18O level.")
    print("Conclusion: To keep the lake's 18O level low, the effect of evaporation must be minimized. This occurs in a 'Wet Climate' with high inflow and 'High Lake Levels', which create a large, stable reservoir.")
    print("-" * 30)

    print("Step 4: Synthesize the Factors")
    print("To form sulfate with low 18O levels, two conditions must be met:")
    print("1. The water entering the lake must be low in 18O. This requires a 'Cold Climate'.")
    print("2. The lake water must not become enriched in 18O via evaporation. This requires a 'Wet Climate' and 'High Lake Levels'.")
    print("\nFinal combination: Wet, cold climate with high lake levels.")
    print("This matches answer choice G.")

# Execute the reasoning and print the final answer in the required format
if __name__ == "__main__":
    solve_paleoclimate_question()
    sys.stdout.write("<<<G>>>\n")